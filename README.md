# **Snakemake workflow: variant calling using the Genome Analysis Toolkit (GATK) best practices**

![GitHub release](https://img.shields.io/github/v/release/kevin-wamae/variant-calling-with-Snakemake-and-GATK)
![GitHub](https://img.shields.io/github/license/kevin-wamae/variant-calling-with-Snakemake-and-GATK?color=maroon)
[![conda](https://img.shields.io/badge/conda->=23.1.0-brightgreen.svg)](https://github.com/conda/conda)
[![snakemake](https://img.shields.io/badge/snakemake-7.24.2-brightgreen.svg)](https://snakemake.readthedocs.io)



[![trimmomatic](https://img.shields.io/badge/trimmomatic-0.39-red.svg)](http://www.usadellab.org/cms/?page=trimmomatic)
[![fastp](https://img.shields.io/badge/fastp-0.23.2-red.svg)](https://github.com/OpenGene/fastp)
[![bwa](https://img.shields.io/badge/bwa-0.7.17-red.svg)](https://github.com/lh3/bwa)
[![samtools](https://img.shields.io/badge/samtools-1.16.1-red.svg)](https://github.com/samtools/samtools)
[![gatk4](https://img.shields.io/badge/gatk4-4.4.0.0-red.svg)](https://github.com/broadinstitute/gatk)
[![samblaster](https://img.shields.io/badge/samblaster-0.1.26-red.svg)](https://github.com/GregoryFaust/samblaster)
[![bcftools](https://img.shields.io/badge/bcftools-1.17-red.svg)](https://github.com/samtools/bcftools)
[![snpeff](https://img.shields.io/badge/snpeff-5.1-red.svg)](http://pcingola.github.io/SnpEff/)


---

## **Table of contents**
- [**Snakemake workflow: variant calling using the Genome Analysis Toolkit (GATK) best practices**](#snakemake-workflow-variant-calling-using-the-genome-analysis-toolkit-gatk-best-practices)
  - [**Table of contents**](#table-of-contents)
  - [**Motivation**](#motivation)
  - [**Pipeline sections**](#pipeline-sections)
      - [  **Step 1 - Compile List of Output Files**](#step-1---compile-list-of-output-files)
      - [  **Step 2 - Gather Genome Data**](#step-2---gather-genome-data)
      - [  **Step 3 - Perform Fastq Quality Control**](#step-3---perform-fastq-quality-control)
      - [  **Step 4 - Map Reads to Genome**](#step-4---map-reads-to-genome)
      - [  **Step 5 - Generate Mapping Quality Statistics**](#step-5---generate-mapping-quality-statistics)
      - [  **Step 6 - Perform Variant Calling**](#step-6---perform-variant-calling)
      - [  **Step 7 - Perform Variant Filtering (Hard or Soft)**](#step-7---perform-variant-filtering-hard-or-soft)
      - [  **Step 8 - Annotatate Variants and Calculate Allele Frequencies**](#step-8---annotatate-variants-and-calculate-allele-frequencies)
  - [**Project dependencies:**](#project-dependencies)
  - [**Where to start**](#where-to-start)
  - [**Directory structure**](#directory-structure)
  - [**Running the analysis**](#running-the-analysis)
  - [**Feedback and Issues**](#feedback-and-issues)


## **Motivation**

- This repository contains a pipeline built with [Snakemake](https://snakemake.readthedocs.io/en/stable/) for variant calling using Illumina-generated sequences and is based on the [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices) for variant calling using either hard or soft filters.
- Additionally, this pipeline aims to reproduce a recently published pipeline that optimized the GATK4 variant calling pipeline for _Plasmodium falciparum_ ([_preprint_](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9949269/)). However, this is not limited to _P. falciparum_ and can be used for any organism of interest.


## **Pipeline sections**
- The pipeline handles paired-end reads and below are the analysis sections in the Snakefile:

#### &nbsp;&nbsp;**Step 1 - Compile List of Output Files**
  - `rule all` - gather all output files

#### &nbsp;&nbsp;**Step 2 - Gather Genome Data**
  - `gather_genome_data`: aggregate genome data from the snpeff folder
  - `gatk_genome_dict`: create genome dictionary for gatk tools
  - `samtools_index`: index the genome fasta file
  - `bedops_gff2bed`: convert the genome annotation .gff to .bed file

#### &nbsp;&nbsp;**Step 3 - Perform Fastq Quality Control**
  -  `trim_reads`: trim adapters and low quality bases using trimmomatic or fastp

#### &nbsp;&nbsp;**Step 4 - Map Reads to Genome**
 - `bwa_index`: generate bwa genome-index files for mapping reads
 - `bwa_mem`: map reads to genome, fixmate, convert .sam to .bam and finally remove artifacts
 - `mark_duplicates`: mark duplicate reads using gatk MarkDuplicatesSpark or Samblaster

#### &nbsp;&nbsp;**Step 5 - Generate Mapping Quality Statistics**
  - `samtools_idxstats`: calculate alignment statistics based on the reference sequence
  - `samtools_flagstats`: calculates and summarizes various alignment statistics
  - `samtools_depth`: calculate the depth of coverage for each position in the genome
  - `gatk_insert_size_metrics`: collect insert size metrics
  - `gatk_alignment_summary_metrics`: generate a summary of alignment metrics from the BAM file

#### &nbsp;&nbsp;**Step 6 - Perform Variant Calling**
  - `gatk_haplotypecaller`: call snps and indels via local re-assembly of haplotypes and generate gVCFs
  - `generate_sample_name_map`: generate a map of sample names and the respective vcf files
  - `gatk_genomics_db_import`: merge gVCFs into one genomic database
  - `gatk_genotype_gvcfs`: perform joint genotyping and generate the final VCF in which all samples have been jointly genotyped

#### &nbsp;&nbsp;**Step 7 - Perform Variant Filtering (Hard or Soft)**
  - `bcftools_normalize`: normalize indels, left-align variants, split multiallelic sites into multiple rows and recover multiallelics from multiple rows
  - hard_filter_variants:
    - `gatk_split_variants`: separate snps and indels into separate vcf files
    - `gatk_filter_hard`: apply hard filters to snps and indels
    - `gatk_merge_vcfs`: merge snps and indels into one vcf file
  - soft_filter_variants:
    - `gatk_vqsr_indels`: perform variant quality score recalibration to indels
    - `gatk_apply_vqsr_indels`: apply variant quality score recalibration to indels
    - `gatk_vqsr_snps`: perform variant quality score recalibration to snps
    - `gatk_apply_vqsr_snps`: apply variant quality score recalibration to snps
  - `gatk_filter_pass`: filter out variants that do not pass the hard or soft filters

#### &nbsp;&nbsp;**Step 8 - Annotatate Variants and Calculate Allele Frequencies**
  - `snpeff_annotate_variants`: variant annotation and functional effect prediction
  - `gatk_variants_to_table`: extract variant information into a table
  - `python` - calculate allele frequencies and transform the summary table from wide to long format

---

## **Project dependencies:**

- [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) - an open-source package management system and environment management system that runs on various platforms, including Windows, MacOS, Linux

- [Snakemake](https://github.com/snakemake/snakemake) - a workflow management system that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment, together with a clean and modern specification language in python style.

---

## **Where to start**

- Install conda for your operating System (_the pipeline is currently tested on Linux and MacOS_):
  - [Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
  - [MacOS](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html)
- Clone this project using the following command in your terminal:
  - `git clone https://github.com/kevin-wamae/variant-calling-with-Snakemake-and-GATK.git`
- Type the following command in your terminal to navigate into the cloned directory using the command below. This will be the root directory of the project:
  - `cd variant-calling-with-Snakemake-and-GATK`
  
- **_Note: All subsequent commands should be run from the root directory of this project. However, users can modify the scripts to their liking_**
 
 ---

## **Directory structure**
- Below is the default directory structure:
    - **config/**   - contains the Snakemake-configuration files
    - **input/** - contains input files
      - **bed/** - contains the bed files for specifying the intervals of interest
      - **fastq/** - contains the FastQ files
      - **known_sites/** - contains the positive-training dataset for variant filtering
    - **output/** - contains numbered-output directories from the analysis
    - **workflow/** - contains the Snakemake workflow files
      - **envs/** - contains the Conda environment-configuration files
      - **scripts/** - contains the scripts used in the pipeline
```
.
|-- config
|-- input
|   |-- bed
|   |-- fastq
|   -- known_sites
|-- output
|-- workflow
    |-- envs
    |-- scripts
```


- This pipeline uses _`global_wildcards()`_ to match FastQ sample names and mate files in the `input/fastq/` directory, using the naming convention below:
  - `reads_R1.fastq.gz` = first mate
  - `reads_R2.fastq.gz` = second mate
  - If you have a different naming convention (eg. [_this_](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm)), you can rename the FastQ files by executing the python script in the _workflow/scripts/_ directory:
    - `python workflow/scripts/fastq_rename.py`
  - Therefore, the user can deposit their FastQ files in the `input/fastq/` directory or edit the `config/config.yaml` file to point to the directory with FastQ files and the pipeline will automatically match the sample names and mates files

- The configuration file (`config/config.yaml`) specifies additional resources and can be modified to suit one's needs, such as:
  - Input files
  - Output directories,
  - The option to choose between tools and methods, e.g.:
    - `fastp` or `trimmomatic` for read trimming
    - `gatk MarkDuplicatesSpark` or `samblaster` for marking duplicates
    - hard or soft filtering of variants
  - Other parameters, such as the number of threads to use

---

## **Running the analysis**
After navigating into the root directory of the project, run the analysis by executing the following commands in your terminal to:

1. Create a conda analysis environment by running the command below in your terminal. This will create a conda environment named `variant-calling-gatk` and install [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [SnpEff](https://pcingola.github.io/SnpEff/se_introduction/) and [Graphviz](https://anaconda.org/anaconda/graphviz) (for visualizing the workflow) in the environment. **_Note:_** This only needs to be done once.
     - `conda env create --file workflow/envs/environment.yaml`

2. Activate the conda environment by running the command below in your terminal. **_Note:_** This needs to be done every time you exit and restart your terminal and want re-run this pipeline
     - `conda activate variant-calling-gatk`

3. Execute the shell script below to create the SnpEff database for variant annotation. This will download the _P. falciparum_ genome data from [PlasmoDB](https://plasmodb.org/) and create a database in the **output/** directory. **_Note:_** This is an important step because the genome-FASTA and GFF files are required for read-mapping and variant calling. It can also be modified to suit one's needs such as download genome files for your organism of interest:
     - `bash workflow/scripts/create_snpeff_db.sh`

4. Finally, execute the whole Snakemake pipeline by running the following command in your terminal:
    - `snakemake --use-conda --cores 2 --jobs 1`
    - This will run the whole pipeline using a maximum of two cores and one job in parallel. The `--cores` flag specifies the number of cores to use for each job and the `--jobs` flag specifies the number of jobs to run in parallel.
    - If you want to run the pipeline using more resources, you can increase the number of cores and jobs. For example, to run the pipeline using 4 cores and 2 jobs in parallel, run the following command:
        - `snakemake --use-conda --cores 4 --jobs 2`
    - Additionally, you can change the `threads` entry in `line 3` of the configuration file (`config/config.yaml`) to specify the number of cores to use for each step in the pipeline.

5. Once the analysis is complete, look through **output/** directory to view the results of the analysis

6. Summary statistics can be generated with stand alone scripts in the `workflow/scripts/` directory:
   - To do this, create an conda environment with the following command:
     - `conda env create --file workflow/envs/variant-calling-stats.yaml`
     - activate the conda environment by running the following command: `conda activate variant-calling-stats`
   - To generate a summary of the raw reads, run the following command and look through the `stats_1_raw_fastq.tsv` file in the project directory:
     - `python workflow/scripts/get_raw_fastq_stats.py`
   - To generate a summary of the trimmed reads, run the following command and look through the `stats_2_trimmed_fastq.tsv` file:
     - `python workflow/scripts/get_trimmed_fastq_stats.py`
   - To generate a summary of the mapped reads, run the following command and look through the `stats_3_mapped_reads.tsv` file:
     - `python workflow/scripts/get_mapped_reads_stats.py`
   - To generate a summary of the variants called, run the following command and look through the `stats_4_variant_calling.tsv` file:
     - `bash workflow/scripts/get_variant_calling_stats.sh`
   - Exit this conda environment by running the following command:
     - `conda deactivate variant-calling-stats`

7. Finally, you can deactivate the variant calling conda environment by running the following command:
     - `conda deactivate variant-calling-gatk`

---

## **Feedback and Issues**

Report any issues or bugs by openning an issue [here](https://github.com/kevin-wamae/variant-calling-with-Snakemake-and-GATK/issues) or contact me via email at **wamaekevin[at]gmail.com**
  
 
