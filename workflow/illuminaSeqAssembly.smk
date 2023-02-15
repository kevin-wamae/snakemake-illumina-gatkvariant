# A Snakemake pipeline for variant calling from illumina sequences
# *********************************************************************


# dependencies
# *********************************************************************
# configuration file
configfile: "workflow/config.yaml"


# global wild cards of sample and matepair list
(SAMPLES,) = glob_wildcards(config["input"]["fastq"] + "{sample}_R1.fastq.gz")

print(SAMPLES)


# all output out
# *********************************************************************
rule all:
    input:
        # ------------------------------------
        # trim_fastq
        expand(config["fastp"]["dir"] + "{sample}_R1.fastq.gz", sample=SAMPLES),
        expand(config["fastp"]["dir"] + "{sample}_R2.fastq.gz", sample=SAMPLES),
        expand(config["fastp"]["dir"] + "{sample}.json", sample=SAMPLES),
        expand(config["fastp"]["dir"] + "{sample}.html", sample=SAMPLES),
        # ------------------------------------
        # get_genome_data
        config["get_genome_data"]["fasta"],
        config["get_genome_data"]["fasta_idx"],
        config["get_genome_data"]["gff"],
        config["get_genome_data"]["regions"],
        # ------------------------------------
        # bwa_index_genome
        config["bwa"]["index"],
        # ------------------------------------
        # bwa_map_reads
        expand(config["bwa"]["dir"] + "{sample}.bam", sample=SAMPLES),
        expand(config["bwa"]["dir"] + "{sample}.bam.bai", sample=SAMPLES),
        expand(config["bwa"]["dir_stats"] + "{sample}.bam.idxstats.txt", sample=SAMPLES),
        expand(
            config["bwa"]["dir_stats"] + "{sample}.bam.flagstats.txt", sample=SAMPLES
        ),
        # bam_mark_dup=config["bwa"]["dir"] + "{sample}.markdup.bam",
        # bam_index=config["bwa"]["dir"] + "{sample}.markdup.bam.bai",
        # bam_to_bed=config["bwa"]["dir"] + "{sample}.markdup.bam.bed",


# genome data - download genome data
# *********************************************************************
rule get_genome_data:
    input:
        genome=config["input"]["genome"]["fasta"],
        gff=config["input"]["genome"]["gff"],
    output:
        genome=config["get_genome_data"]["fasta"],
        genome_index=config["get_genome_data"]["fasta_idx"],
        gff=config["get_genome_data"]["gff"],
        regions=config["get_genome_data"]["regions"],
    params:
        loci=config["get_genome_data"]["loci"],
        feature=config["get_genome_data"]["feature"],
    run:
        shell(  # copy genome fasta file from snpeff database location
            """
            cp -f {input.genome} {output.genome}
            """
        )
        shell(  # faidx - generate fasta index file
            """
            samtools faidx {output.genome}
            """
        )
        shell(  # copy annotation file from snpeff database location
            """
            cp -f {input.gff} {output.gff}
            """
        )
        shell(  # extract regions of interest, here we extract AMA1 and K13 coding regions
            """
            grep {output.gff} -e '{params.loci}' |\
            grep -e '{params.feature}' |\
            convert2bed --input=gff --output=bed > {output.regions}
            """
        )


# fastp - clip illumina adapters, paired end mode
# *********************************************************************
rule trim_fastq_files:
    input:
        in1=config["input"]["fastq"] + "{sample}_R1.fastq.gz",
        in2=config["input"]["fastq"] + "{sample}_R2.fastq.gz",
    output:
        out1=config["fastp"]["dir"] + "{sample}_R1.fastq.gz",
        out2=config["fastp"]["dir"] + "{sample}_R2.fastq.gz",
        json=config["fastp"]["dir"] + "{sample}.json",
        html=config["fastp"]["dir"] + "{sample}.html",
    params:
        threads=config["extra"]["threads"],
        min_len=config["fastp"]["min_len"],
    shell:
        """
        fastp \
            --thread {params.threads} \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --cut_tail --cut_tail_mean_quality 20 \
            --length_required {params.min_len} \
            --in1 {input.in1} \
            --in2 {input.in2} \
            --out1 {output.out1} --out2 {output.out2} \
            --json {output.json} --html {output.html}
            """


# bwa - generate bwa genome-index files ------
# *********************************************************************
rule bwa_index_genome:
    input:
        genome=rules.get_genome_data.output.genome,
    output:
        genomeIndex=touch(config["bwa"]["index"]),
    shell:
        """
        bwa index -p {output.genomeIndex} {input.genome}
        """


# bwa/samtools/sambamba: [-a bwtsw|is]
# *********************************************************************
rule bwa_map_reads:
    input:
        genome=rules.get_genome_data.output.genome,
        genomeIndex=rules.bwa_index_genome.output.genomeIndex,
        fastqR1=rules.trim_fastq_files.output.out1,
        fastqR2=rules.trim_fastq_files.output.out1,
        regions=rules.get_genome_data.output.regions,
    output:
        bam=config["bwa"]["dir"] + "{sample}.bam",
        index=config["bwa"]["dir"] + "{sample}.bam.bai",
        idxstats=config["bwa"]["dir_stats"] + "{sample}.bam.idxstats.txt",
        flagstats=config["bwa"]["dir_stats"] + "{sample}.bam.flagstats.txt",
        # bam_to_bed=config["bwa"]["dir"] + "{sample}.markdup.bam.bed",
    params:
        threads=config["extra"]["threads"],
    run:
        shell(  # bwa mem - map reads to reference genome
            """
            bwa mem -M -t \
                {params.threads} \
                {input.genomeIndex} \
                {input.fastqR1} {input.fastqR2} |\
            samblaster |\
            samtools view -S -b  --targets-file {input.regions} |\
            samtools sort -o {output.bam}
            """
        )
        shell(  # samtools index - index bam file for downstream analysis
            """
            samtools index {output.bam} {output.index}
            """
        )
        shell(  # stats for the bam index file
            """
            samtools idxstats {output.bam} > {output.idxstats}
            """
        )
        shell(  # counts for each of 13 categories based primarily on bit flags in the FLAG field
            """
            samtools flagstats {output.bam} --output-fmt tsv > {output.flagstats}
            """
        )


# # *********************************************************************
# # bedops - convert .gff to .bed using bedops and filter to exome
# rule gff_to_bed:
#     input:
#         gff=config["input"]["genome"]["gff"],
#     output:
#         bed=config["bed"]["gffToBed"],
#     params:
#         gffFilter="protein_coding_gene",
#     shell:
#         """
#        convert2bed --input=gff --output=bed < {input.gff} | \
#             grep -e {params.gffFilter} > {output.bed}
#         """
# # *********************************************************************
# # bedtools - merge overlapping intervals
# #  -s: merge features that are on the same strand
# #  -i: input (bed/gff/vcf)
# #  -c: report strand from bed file
# #  -o: specify operation (distinct removes duplicates)
# rule merge_features:
#     input:
#         bed=rules.map_reads.output.bam_to_bed,
#         bam=rules.map_reads.output.bam_mark_dup,
#     output:
#         merged=config["bed"]["dir"] + "{sample}.merged.bed",
#         bedCov=config["bed"]["dir"] + "{sample}.coverage.bed",
#         annotated=config["bed"]["dir"] + "{sample}.annotated.bed",
#     params:
#         bed=rules.gff_to_bed.output.bed,
#     run:
#         shell(  # bedtools - merge overlapping intervals
#             """
#             bedtools merge -s -c 6,5 -o distinct,mean -i {input.bed} > {output.merged}
#             """
#         )
#         shell(  # bedtools - compute coverate across genome features
#             """
#             bedtools genomecov -bg -ibam {input.bam} |\
#             bedtools merge -c 4 -o median > {output.bedCov}
#             """
#         )
#         shell(  # bedtools - annotate intervals in the bed file with features from
#             # column 10 of .gff and collapse overlapping features
#             """
#             bedtools map -c 10 -o collapse -a {output.bedCov} -b {params.bed} > {output.annotated}
#             """
#         )
# # *********************************************************************
# # bcftools - variant calling
# #   - bcftools mpileup - generates genotype likelihoods at each genomic
# #     position with coverage.
# #   - bcftools call - makes the actual calls
# #   - bcftools filter - drop variants with QUAL<=20 and Depth of Coverage
# rule bcftools_variant_calling:
#     input:
#         genomeIndex=rules.create_genome_index.input.genomeFasta,
#         bam=rules.map_reads.output.bam_mark_dup,
#         regions=config["input"]["genome"]["regions"],
#     output:
#         vcf=config["bcftools"]["dir"] + "{sample}.vcf.gz",
#     params:
#         threads=config["extra"]["threads"],
#     shell:
#         """
#         bcftools mpileup \
#             --threads {params.threads} \
#             --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,FORMAT/QS,INFO/AD,INFO/ADF,INFO/ADR \
#             --fasta-ref {input.genomeIndex} {input.bam} \
#             --output-type z |\
#         bcftools call \
#             --threads {params.threads} \
#             --regions {input.regions} \
#             --skip-variants indels \
#             --multiallelic-caller --variants-only |\
#         bcftools filter \
#             --threads {params.threads} \
#             --soft-filter LowQual \
#             --include 'QUAL>20 || DP>100' \
#             --output-type z --output {output.vcf}
#         """
# # *********************************************************************
# # snpEff - variant annotation and functional effect prediction
# # *********************************************************************
# rule snpeff_annotate_variants:
#     input:
#         rules.bcftools_variant_calling.output.vcf,
#     output:
#         vcf=config["snpEff"]["dir"] + "{sample}.vcf.gz",
#     params:
#         config=config["snpEff"]["config"],
#         database=config["snpEff"]["database"],
#     shell:
#         """
#         snpEff -no-downstream -no-intergenic -no-intron -no-upstream \
#                -no-utr -hgvs1LetterAa -noLof -noShiftHgvs -noMotif \
#                -no SPLICE_SITE_REGION -noInteraction -noStats  \
#                -config {params.config} {params.database}  {input} | gzip > {output.vcf}
#         """
# # *********************************************************************
# # SnpSift - extract vcf fields
# # *********************************************************************
# rule snpsift_extract_variants:
#     input:
#         rules.snpeff_annotate_variants.output.vcf,
#     output:
#         allele_list=config["snpSift"]["dir"] + "{sample}.allele.txt",
#         allele_freq=config["snpSift"]["dir"] + "{sample}.allele.freq.txt",
#     run:
#         shell(  # SnpSift - extract vcf fields
#             """
#                  SnpSift extractFields {input} \
#                     CHROM POS REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" \
#                     "ANN[*].GENEID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" \
#                     "ANN[*].CDS_POS" "ANN[*].AA_POS" "GEN[*].AD" > {output.allele_list}
#                 """
#         )
#         shell(  # awk - compute allele frequencies from AD column of .vcf file
#             """
#             awk '
#                 BEGIN {{ FS=OFS="\\t" }}
#                 NR == 1 {{
#                     allelFreq1 = "AF_REF"
#                     allelFreq2 = "AF_ALT"
#                 }}
#                 NR > 1 {{
#                     split($12,a,",")
#                     sum = a[1] + a[2]
#                     if ( sum ) {{
#                         allelFreq1 = a[1] / sum
#                         allelFreq2 = a[2] / sum
#                     }}
#                     else {{
#                         allelFreq1 = 0
#                         allelFreq2 = 0
#                     }}
#                 }}
#                 {{ print $0, allelFreq1, allelFreq2 }}
#             ' {output.allele_list} | sed -e 's/ANN\\[\\*\\]\\.\\|GEN\\[\\*\\]\\.//g' > {output.allele_freq}
#             """
#         )
