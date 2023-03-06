# #####################################################################
#   A Snakemake pipeline for variant calling from illumina sequences
# #####################################################################


# dependencies
# *********************************************************************
# configuration file
configfile: "workflow/config.yaml"


# global wild cards of sample and matepair list
# - this assumes that all fastq files have the SampleName_R1.fastq.gz and SampleName_R2.fastq.gz format
# - If not, they might have the SampleName_S1_L001_R1_001.fastq.gz and SampleName_S1_L001_R2_001.fastq.gz format
# - In this case, rename the files running the following command in your terminal, at the top level of the project:
# python3 workflow/rename_fastq_files.py
(SAMPLES,) = glob_wildcards(config["input"]["fastq"] + "{sample}_R1.fastq.gz")


# a list of all output files
# *********************************************************************
rule all:
    input:
        # ------------------------------------
        # get_genome_data
        config["get_genome_data"]["fasta"],
        config["get_genome_data"]["gff"],
        # config["get_genome_data"]["regions"],
        # ------------------------------------
        # samtools_index
        config["samtools_index"]["fasta_idx"],
        # ------------------------------------        
        # bedops_gff2bed
        config["bedops_gff2bed"]["bed"],
        # # ------------------------------------
        # # trim_fastq
        # expand(config["fastp"]["dir"] + "{sample}_R1.fastq.gz", sample=SAMPLES),
        # expand(config["fastp"]["dir"] + "{sample}_R2.fastq.gz", sample=SAMPLES),
        # expand(config["fastp"]["log"] + "{sample}.json", sample=SAMPLES),
        # expand(config["fastp"]["log"] + "{sample}.html", sample=SAMPLES),
        # # ------------------------------------
        # # bwa_index_genome
        # config["bwa"]["genome_index"],
        # # ------------------------------------
        # # bwa_map_reads
        # expand(config["bwa"]["dir"] + "{sample}.bam", sample=SAMPLES),
        # expand(config["bwa"]["dir"] + "{sample}.bam.bai", sample=SAMPLES),
        # # ------------------------------------
        # # samtools_mapping_stats
        # expand(
        #     config["mapping_stats"]["dir"] + "{sample}.bam.idxstats.txt",
        #     sample=SAMPLES,
        # ),
        # expand(
        #     config["mapping_stats"]["dir"] + "{sample}.bam.flagstats.txt",
        #     sample=SAMPLES,
        # ),
        # # ------------------------------------
        # # bcftools_variant_calling
        # expand(config["bcftools"]["dir"] + "{sample}.vcf.gz", sample=SAMPLES),
        # # ------------------------------------
        # # snpeff_annotate_vcf
        # expand(config["snpEff"]["dir"] + "{sample}.vcf.gz", sample=SAMPLES),
        # # ------------------------------------
        # # snpsift_filter_vcf
        # expand(config["snpSift"]["dir"] + "{sample}.allele.txt", sample=SAMPLES),
        # expand(config["snpSift"]["dir"] + "{sample}.allele.freq.txt", sample=SAMPLES),


# genome data - download genome data
# *********************************************************************
rule get_genome_data:
    input:
        genome=config["input"]["genome"]["fasta"],
        gff=config["input"]["genome"]["gff"],
    output:
        genome=config["get_genome_data"]["fasta"],
        gff=config["get_genome_data"]["gff"],
        # bed=config["get_genome_data"]["bed"],
        # regions=config["get_genome_data"]["regions"],
    params:
        # loci=config["get_genome_data"]["loci"],
        # feature_genome=config["get_genome_data"]["feature_filter_genome"],
        # feature_variants=config["get_genome_data"]["feature_filter_variants"],
    run:
        shell(  # cp - copy genome fasta file from snpeff database location
            """
            cp -f {input.genome} {output.genome}
            """
        )
        shell(  # cp - copy annotation file from snpeff database location
            """
            cp -f {input.gff} {output.gff}
            """
        )
        # shell(
        #     # bedops - convert genome annotation GFF to BED
        #     #   - compared to the BED file containing the regions of interest (below), we will select
        #     #     entire protein coding regions of the genome
        #     """
        #     convert2bed --input=gff --output=bed < {input.gff} | \
        #         grep -e {params.feature_genome} > {output.bed}
        #     """
        # )
        # shell(
        #     # grep and bedops - extract regions of interest, here we extract AMA1 and K13 coding regions
        #     #   - compared to the BED file containing the protein coding regions (above), we will only
        #     #     select regions within the coding regions of the genome
        #     """
        #     grep {output.gff} -e '{params.loci}' |\
        #     grep -e '{params.feature_variants}' |\
        #     convert2bed --input=gff --output=bed > {output.regions}
        #     """
        # )


# samtools index - index genome fasta file
# *********************************************************************s
rule samtools_index:
    input:
        rules.get_genome_data.output.genome,
    output:
        config["samtools_index"]["fasta_idx"],
    wrapper:
        "master/bio/samtools/faidx"


# bedops - convert genome annotation GFF to BED
# *********************************************************************
rule bedops_gff2bed:
    input:
        rules.get_genome_data.output.gff,
    output:
        config["bedops_gff2bed"]["bed"],
    params:
        feature=config["bedops_gff2bed"]["feature"],
    conda:
        config["conda_env"]["bedops"],
    shell:
        """
        convert2bed --input=gff --output=bed < {input} | \
            grep -e {params.feature} > {output}
        """


# # fastp - clip illumina adapters, paired end mode
# # *********************************************************************
# rule trim_fastq_files:
#     input:
#         in1=config["input"]["fastq"] + "{sample}_R1.fastq.gz",
#         in2=config["input"]["fastq"] + "{sample}_R2.fastq.gz",
#     output:
#         out1=config["fastp"]["dir"] + "{sample}_R1.fastq.gz",
#         out2=config["fastp"]["dir"] + "{sample}_R2.fastq.gz",
#         json=config["fastp"]["log"] + "{sample}.json",
#         html=config["fastp"]["log"] + "{sample}.html",
#     params:
#         threads=config["extra"]["threads"],
#         min_len=config["fastp"]["min_len"],
#         min_qual=config["fastp"]["min_qual"],
#     shell:
#         """
#         fastp \
#             --thread {params.threads} \
#             --detect_adapter_for_pe \
#             --length_required {params.min_len} \
#             --qualified_quality_phred 20 \
#             --cut_tail --cut_tail_mean_quality {params.min_qual} \
#             --in1 {input.in1} \
#             --in2 {input.in2} \
#             --out1 {output.out1} --out2 {output.out2} \
#             --json {output.json} --html {output.html}
#             """


# # #####################################################################
# #                      BCFTOOLS VARIANT CALLING
# # # ###################################################################


# # bwa - generate bwa genome-index files
# # *********************************************************************
# rule bwa_index_genome:
#     input:
#         genome=rules.get_genome_data.output.genome,
#     output:
#         genome_index=touch(config["bwa"]["genome_index"]),
#     shell:
#         """
#         bwa index -p {output.genome_index} {input.genome}
#         """


# # bwa/samtools/sambamba: [-a bwtsw|is]
# # *********************************************************************
# rule bwa_map_reads:
#     input:
#         genome=rules.get_genome_data.output.genome,
#         genome_index=rules.bwa_index_genome.output.genome_index,
#         fastqR1=rules.trim_fastq_files.output.out1,
#         fastqR2=rules.trim_fastq_files.output.out2,
#         regions=rules.get_genome_data.output.regions,
#     output:
#         bam=config["bwa"]["dir"] + "{sample}.bam",
#         index=config["bwa"]["dir"] + "{sample}.bam.bai",
#     params:
#         threads=config["extra"]["threads"],
#         mapping_qual=config["bwa"]["mapping_qual"],
#         exclude_flag=config["bwa"]["exclude_flag"],
#         rg_id="{sample}",
#         rg_sm="{sample}",
#         rg_lb=config["bwa"]["read_group"]["library"],
#         rg_pl=config["bwa"]["read_group"]["platform"],
#     run:
#         shell(
#             # bwa mem - map reads to reference genome
#             # samblaster - mark duplicates
#             # samtools view - convert sam to bam
#             #  - remove reads with mapping quality < 60
#             #  - remove reads with flag 4 (unmapped)
#             # samtools sort - sort bam file
#             """
#             bwa mem \
#                 -M \
#                 -t {params.threads} \
#                 -R '@RG\\tID:{params.rg_id}\\tSM:{params.rg_sm}\\tLB:{params.rg_lb}\\tPL:{params.rg_pl}' \
#                 {input.genome_index} \
#                 {input.fastqR1} {input.fastqR2} |\
#             samblaster -M \
#                 --removeDups \
#                 --addMateTags |\
#             samtools view -S -b \
#                 --targets-file {input.regions} \
#                 --min-MQ {params.mapping_qual} \
#                 --exclude-flags {params.exclude_flag} |\
#             samtools sort -o {output.bam}
#             """
#         )
#         shell(  # samtools index - index bam file for downstream analysis
#             """
#             samtools index {output.bam} {output.index}
#             """
#         )
# # get mapping-quality statistics from BAM file
# # *********************************************************************
# rule samtools_mapping_stats:
#     input:
#         bam=rules.bwa_map_reads.output.bam,
#     output:
#         idxstats=config["mapping_stats"]["dir"] + "{sample}.bam.idxstats.txt",
#         flagstats=config["mapping_stats"]["dir"] + "{sample}.bam.flagstats.txt",
#     run:
#         shell(  # samtools idxstats - counts for each of 13 categories based primarily on bit flags in the FLAG field
#             """
#             samtools idxstats {input.bam} > {output.idxstats}
#             """
#         )
#         shell(  # samtools flagstats - stats for the bam index file
#             """
#             samtools flagstats {input.bam} --output-fmt tsv > {output.flagstats}
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
#         genome=rules.get_genome_data.output.genome,
#         bam=rules.bwa_map_reads.output.bam,
#     output:
#         vcf=config["bcftools"]["dir"] + "{sample}.vcf.gz",
#     params:
#         threads=config["extra"]["threads"],
#     shell:
#         """
#         bcftools mpileup \
#             --threads {params.threads} \
#             --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,FORMAT/QS,INFO/AD,INFO/ADF,INFO/ADR \
#             --fasta-ref {input.genome} \
#             {input.bam} |\
#         bcftools call \
#             --threads {params.threads} \
#             --skip-variants indels \
#             --multiallelic-caller \
#             --variants-only |\
#         bcftools filter \
#             --threads {params.threads} \
#             --exclude 'QUAL<20' |\
#         bcftools view \
#             --threads {params.threads} \
#             --include 'DP>20' \
#             --output-type z --output {output.vcf}
#         """
# # *********************************************************************
# # snpEff - variant annotation and functional effect prediction
# # *********************************************************************
# rule snpeff_annotate_vcf:
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
#             -no-utr -hgvs1LetterAa -noLof -noShiftHgvs -noMotif \
#             -no SPLICE_SITE_REGION -noInteraction -noStats  \
#             -config {params.config} {params.database} {input} | gzip > {output.vcf}
#         """
# # *********************************************************************
# # SnpSift - extract vcf fields
# # *********************************************************************
# rule snpsift_filter_vcf:
#     input:
#         rules.snpeff_annotate_vcf.output.vcf,
#     output:
#         allele_list=config["snpSift"]["dir"] + "{sample}.allele.txt",
#         allele_freq=config["snpSift"]["dir"] + "{sample}.allele.freq.txt",
#     run:
#         shell(  # SnpSift - extract vcf fields
#             """
#             SnpSift extractFields {input} \
#                 CHROM POS REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" \
#                 "ANN[*].GENEID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" \
#                 "ANN[*].CDS_POS" "ANN[*].AA_POS" "GEN[*].GT" "GEN[*].AD" > {output.allele_list}
#             """
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