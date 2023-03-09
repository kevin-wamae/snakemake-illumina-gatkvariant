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
        # gather_genome_data
        config["gather_genome_data"]["fasta"],
        config["gather_genome_data"]["gff"],
        # ------------------------------------
        # genome_dict
        config["gather_genome_data"]["dict"],
        # config["gather_genome_data"]["regions"],
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
        # ------------------------------------
        # gatk_markdup
        expand(config["gatk_markdup"]["dir"] + "{sample}.bam", sample=SAMPLES),
        expand(
            config["gatk_markdup"]["metrics"] + "{sample}.metrics.txt", sample=SAMPLES
        ),
        # ------------------------------------
        # samtools
        expand(config["samtools_view"]["dir"] + "{sample}.bam", sample=SAMPLES),
        expand(config["samtools_view"]["dir"] + "{sample}.bam.bai", sample=SAMPLES),
        # ------------------------------------
        # samtools_idxstats & samtools_flagstats
        expand(
            config["samtools_stats"]["dir"] + "{sample}.bam.idxstats.txt",
            sample=SAMPLES,
        ),
        expand(
            config["samtools_stats"]["dir"] + "{sample}.bam.flagstat.txt",
            sample=SAMPLES,
        ),
        # ------------------------------------
        # gatk_collect_insert_size_metrics
        expand(
            config["gatk_insert_size"]["dir_metrics"] + "{sample}.metrics.txt",
            sample=SAMPLES,
        ),
        expand(
            config["gatk_insert_size"]["dir_histogram"] + "{sample}.histogram.pdf",
            sample=SAMPLES,
        ),
        # ------------------------------------
        # gatk_realign
        expand(
            config["gatk_haplotypecaller"]["dir"] + "{sample}.vcf.gz", sample=SAMPLES
        ),
        # ------------------------------------
        # gatk_haplotypecaller
        expand(
            config["gatk_haplotypecaller"]["dir"] + "{sample}.vcf.gz", sample=SAMPLES
        ),
        # # ------------------------------------
        # gatk_genomics_db_import
        config["gatk_genomicsdbimport"]["dir"],
        # config["gatk_genomicsdbimport"]["dir"],
        # # ------------------------------------
        # # snpeff_annotate_vcf
        # expand(config["snpEff"]["dir"] + "{sample}.vcf.gz", sample=SAMPLES),
        # # ------------------------------------
        # # snpsift_filter_vcf
        # expand(config["snpSift"]["dir"] + "{sample}.allele.txt", sample=SAMPLES),
        # expand(config["snpSift"]["dir"] + "{sample}.allele.freq.txt", sample=SAMPLES),


# ######################################################################
#                     step 1 - prepare genome data
# ######################################################################


# genome data - download genome data
# *********************************************************************
rule gather_genome_data:
    input:
        genome=config["input"]["genome"]["fasta"],
        gff=config["input"]["genome"]["gff"],
    output:
        genome=config["gather_genome_data"]["fasta"],
        gff=config["gather_genome_data"]["gff"],
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
        genome=rules.gather_genome_data.output.genome,
    output:
        genome_dict=config["gather_genome_data"]["dict"],
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
        rules.gather_genome_data.output.genome,
    output:
        config["samtools_index"]["fasta_idx"],
    wrapper:
        "master/bio/samtools/faidx"


# bedops - convert genome GFF to BED
# *********************************************************************
rule bedops_gff2bed:
    input:
        rules.gather_genome_data.output.gff,
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


# ######################################################################
#                     step 2 - quality control
# ######################################################################


# trimmomatic - clip illumina adapters, paired end mode
# TODO: 1. add support for single end mode
# TODO: 2. review the parameters
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


# ######################################################################
#                            step 3 - mapping
# ######################################################################


# bwa - generate bwa genome-index files for mapping
# *********************************************************************
rule bwa_index:
    input:
        genome=rules.gather_genome_data.output.genome,
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
        genome=rules.gather_genome_data.output.genome,
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
        genome=rules.gather_genome_data.output.genome,
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


# gatk - mark duplicates
# *********************************************************************
rule gatk_markdup:
    input:
        bam=rules.gatk_sort_sam.output.sorted,
    output:
        bam=config["gatk_markdup"]["dir"] + "{sample}.bam",
        metrics=config["gatk_markdup"]["metrics"] + "{sample}.metrics.txt",
    log:
        config["gatk_markdup"]["log"] + "{sample}.log",
    params:
        extra="",
        java_opts="",
        #spark_runner="",  # optional, local by default
        #spark_master="",  # optional
        #spark_extra="", # optional
    threads: 8
    wrapper:
        "master/bio/gatk/markduplicatesspark"


# samtools - view (keep reads in core genome regions of BED file)
# *********************************************************************
rule samtools_view:
    input:
        bam=rules.gatk_markdup.output.bam,
        genome=rules.gather_genome_data.output.genome,
        core_genome=config["gather_genome_data"]["core"],
    output:
        bam=config["samtools_view"]["dir"] + "{sample}.bam",
        index=config["samtools_view"]["dir"] + "{sample}.bam.bai",
    log:
        config["samtools_view"]["log"] + "{sample}.log",
    threads: config["threads"]
    conda:
        config["conda_env"]["samtools"]
    shell:
        """
        samtools view \
            -b \
            -h \
            -@ {threads} \
            -T {input.genome} \
            -L {input.core_genome} \
            {input.bam} \
            > {output.bam}
        samtools index {output.bam} {output.index}
        """


# samtools - idxstats (get mapping-quality statistics from BAM file)
# *********************************************************************
rule samtools_idxstats:
    input:
        rules.samtools_view.output.bam,
    output:
        config["samtools_stats"]["dir"] + "{sample}.bam.idxstats.txt",
    conda:
        config["conda_env"]["samtools"]
    shell:
        """
        samtools idxstats {input} > {output}
        """


# samtools - flagstats (get mapping-quality statistics from BAM file)
# *********************************************************************
rule samtools_flagstat:
    input:
        rules.samtools_view.output.bam,
    output:
        config["samtools_stats"]["dir"] + "{sample}.bam.flagstat.txt",
    conda:
        config["conda_env"]["samtools"]
    shell:
        """
        samtools flagstat {input} > {output}
        """


# gatk CollectInsertSizeMetrics
# *********************************************************************
rule gatk_collect_insert_size_metrics:
    input:
        bam=rules.samtools_view.output.bam,
        genome=rules.gather_genome_data.output.genome,
    output:
        metrics=config["gatk_insert_size"]["dir_metrics"] + "{sample}.metrics.txt",
        histogram=config["gatk_insert_size"]["dir_histogram"] + "{sample}.histogram.pdf",
    log:
        config["gatk_insert_size"]["log"] + "{sample}.log",
    params:
        extra="",
        java_opts="",
    conda:
        config["conda_env"]["gatk"]
    shell:
        """
        gatk CollectInsertSizeMetrics \
            {params.extra} \
            -R {input.genome} \
            -I {input.bam} \
            -O {output.metrics} \
            -H {output.histogram}
        """


# ######################################################################
#                        step 4 - variant calling
# ######################################################################


# gatk HaplotypeCaller - generate gVCFs
# *********************************************************************
rule gatk_haplotypecaller:
    input:
        bam=rules.samtools_view.output.bam,
        genome=rules.gather_genome_data.output.genome,
        intervals=config["gather_genome_data"]["core"],
    output:
        vcf=config["gatk_haplotypecaller"]["dir"] + "{sample}.vcf.gz",
    log:
        config["gatk_haplotypecaller"]["log"] + "{sample}.log",
    params:
        extra=config["gatk_haplotypecaller"]["extra"],
        java_opts="",
        threads=config["threads"],
    threads: config["threads"]
    conda:
        config["conda_env"]["gatk"]
    shell:
        """
        gatk HaplotypeCaller \
            {params.extra} \
            --native-pair-hmm-threads {params.threads} \
            -R {input.genome} \
            -L {input.intervals} \
            -I {input.bam} \
            -O {output.vcf}
        """


# gatk GenomicsDBImport - merge gVCFs into one genomic database
# *********************************************************************
rule gatk_genomics_db_import:
    input:
        vcf=expand(rules.gatk_haplotypecaller.output.vcf, sample=SAMPLES),
        genome=rules.gather_genome_data.output.genome,
        intervals=config["gather_genome_data"]["core"],
    output:
        dir=directory(config["gatk_genomicsdbimport"]["dir"]),
    log:
        config["gatk_genomicsdbimport"]["dir"] + "genomicsdb.log",
    params:
        java_opts="",
        extra=config["gatk_genomicsdbimport"]["extra"],
        threads=config["threads"],
    threads: config["threads"]
    conda:
        config["conda_env"]["gatk"]
    shell:
        """
        gatk GenomicsDBImport \
            {params.extra} \
            --reader-threads {params.threads} \
            -V {input.vcf} \
            --genomicsdb-workspace-path {output.dir} \
            -L {input.intervals}
        """
