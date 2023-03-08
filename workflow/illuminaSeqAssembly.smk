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
        # ------------------------------------
        # genome_dict
        config["get_genome_data"]["dict"],
        # config["get_genome_data"]["regions"],
        # ------------------------------------
        # samtools_index
        config["samtools_index"]["fasta_idx"],
        # ------------------------------------        
        # bedops_gff2bed
        config["bedops_gff2bed"]["bed"],
        # ------------------------------------
        # trimmomaticS,
        expand(config["trimmomatic"]["dir"] + "{sample}_R1.fastq.gz", sample=SAMPLES),
        expand(config["trimmomatic"]["dir"] + "{sample}_R2.fastq.gz", sample=SAMPLES),
        expand(
            config["trimmomatic"]["dir"] + "{sample}_R1.unpaired.fastq.gz",
            sample=SAMPLES,
        ),
        expand(
            config["trimmomatic"]["dir"] + "{sample}_R2.unpaired.fastq.gz",
            sample=SAMPLES,
        ),
        # ------------------------------------
        # bwa_index
        config["bwa"]["index"],
        # ------------------------------------
        # bwa_mem
        expand(config["bwa"]["dir"] + "{sample}.bam", sample=SAMPLES),
        # ------------------------------------
        # gatk_clean
        expand(config["gatk_clean"]["dir"] + "{sample}.bam", sample=SAMPLES),
        # ------------------------------------
        # gatk_sort
        expand(config["gatk_sort"]["dir"] + "{sample}.bam", sample=SAMPLES),
        # # ------------------------------------
        # # gatk_markdup
        # expand(config["gatk_markdup"]["dir"] + "{sample}.bam", sample=SAMPLES),
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


# genome data - download genome data
# *********************************************************************
rule gatk_genome_dict:
    input:
        genome=rules.get_genome_data.output.genome,
    output:
        genome_dict=config["get_genome_data"]["dict"],
    conda:
        config["conda_env"]["gatk"]
    shell:
        """
        gatk CreateSequenceDictionary \
        --REFERENCE {input.genome} \
        --OUTPUT {output.genome_dict}
        """


# samtools index - index genome fasta file
# *********************************************************************s
rule samtools_index:
    input:
        rules.get_genome_data.output.genome,
    output:
        config["samtools_index"]["fasta_idx"],
    wrapper:
        "master/bio/samtools/faidx"


# bedops - convert genome GFF to BED
# *********************************************************************
rule bedops_gff2bed:
    input:
        rules.get_genome_data.output.gff,
    output:
        config["bedops_gff2bed"]["bed"],
    params:
        feature=config["bedops_gff2bed"]["feature"],
    conda:
        config["conda_env"]["bedops"]
    shell:
        """
        convert2bed --input=gff --output=bed < {input} | \
            grep -e {params.feature} > {output}
        """


# trimmomatic - clip illumina adapters, paired end mode
# *********************************************************************
rule trimmomatic:
    input:
        r1=config["input"]["fastq"] + "{sample}_R1.fastq.gz",
        r2=config["input"]["fastq"] + "{sample}_R2.fastq.gz",
    output:
        r1=config["trimmomatic"]["dir"] + "{sample}_R1.fastq.gz",
        r2=config["trimmomatic"]["dir"] + "{sample}_R2.fastq.gz",
        r1_unpaired=config["trimmomatic"]["dir"] + "{sample}_R1.unpaired.fastq.gz",
        r2_unpaired=config["trimmomatic"]["dir"] + "{sample}_R2.unpaired.fastq.gz",
    params:
        trimmer=config["trimmomatic"]["trimmer"],
        extra=config["trimmomatic"]["extra"],
    log:
        config["trimmomatic"]["dir"] + "log/{sample}.log",
    wrapper:
        "master/bio/trimmomatic/pe"


# bwa - generate bwa genome-index files for mapping
# *********************************************************************
rule bwa_index:
    input:
        genome=rules.get_genome_data.output.genome,
    output:
        index=touch(config["bwa"]["index"]),
    conda:
        config["conda_env"]["bwa"]
    shell:
        """
        bwa index -p {output.index} {input.genome}
        """


# bwa - map reads to genome
# *********************************************************************
rule bwa_mem:
    input:
        reads=[
            rules.trimmomatic.output.r1,
            rules.trimmomatic.output.r2,
        ],
        idx=rules.bwa_index.output,
    output:
        config["bwa"]["dir"] + "{sample}.bam",
    log:
        config["bwa"]["log"] + "{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="none",
        sort_extra="",  # Extra args for samtools/picard.
    threads: config["threads"]
    wrapper:
        "master/bio/bwa/mem"


# gatk - clean sam file (remove artifacts in SAM/BAM files)
# *********************************************************************
rule gatk_clean_sam:
    input:
        bam=rules.bwa_mem.output,
        genome=rules.get_genome_data.output.genome,
    output:
        clean=config["gatk_clean"]["dir"] + "{sample}.bam",
    params:
        java_opts="",
    conda:
        config["conda_env"]["gatk"]
    shell:
        """
        gatk CleanSam \
            -R {input.genome} \
            -I {input.bam} \
            -O {output.clean}
        """


# gatk - sort sam
# *********************************************************************
rule gatk_sort_sam:
    input:
        bam=rules.gatk_clean_sam.output.clean,
        genome=rules.get_genome_data.output.genome,
    output:
        sorted=config["gatk_sort"]["dir"] + "{sample}.bam",
    params:
        java_opts="",
    conda:
        config["conda_env"]["gatk"]
    shell:
        """
        gatk SortSam \
            -R {input.genome} \
            -I {input.bam} \
            -O {output.sorted} \
            --SORT_ORDER coordinate
        """


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
