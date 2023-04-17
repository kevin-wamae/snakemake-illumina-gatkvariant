# **Snakemake workflow: variant calling using GATK4 best practices**

## **Table of contents**
- [**Snakemake workflow: variant calling using GATK best practices**](#snakemake-workflow-variant-calling-using-gatk-best-practices)
  - [**Table of contents**](#table-of-contents)
  - [**Motivation**](#motivation)
  - [**Pipeline sections**](#pipeline-sections)
      - [  **Step 1 - output files**](#step-1---output-files)
      - [  **Step 2 - gather genome data**](#step-2---gather-genome-data)
      - [  **Step 3 - quality control**](#step-3---quality-control)
      - [  **Step 4 - map reads to genome**](#step-4---map-reads-to-genome)
      - [  **Step 5 - mapping quality statistics**](#step-5---mapping-quality-statistics)
      - [  **Step 6 - variant calling**](#step-6---variant-calling)
      - [  **Step 7 - variant filtering**](#step-7---variant-filtering)
      - [  **Step 8 - variant annotation and extraction**](#step-8---variant-annotation-and-extraction)
  - [**Project dependencies:**](#project-dependencies)
  - [**Where to start**](#where-to-start)
  - [**Directory structure**](#directory-structure)
  - [**Running the analysis**](#running-the-analysis)
  - [**Feedback and Issues**](#feedback-and-issues)


## **Motivation**

- This repository contains a pipeline built with [Snakemake](https://snakemake.readthedocs.io/en/stable/) for variant calling using Illumina-generated sequences and is based on the [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices) for variant calling.
- Additionally, this pipeline aims to reproduce a recently published pipeline that optimized the GATK4 variant calling pipeline for _Plasmodium falciparum_ ([_preprint_](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9949269/)). However, this is not limited to _P. falciparum_ and can be used for any organism of interest.
- _**Note: The pipeline implements VCF hard-filtering, instead of the recommended soft-filtering via Variant Quality Score Recalibration (VQSR), which will be implemented in a future release**_.


## **Pipeline sections**
- The pipeline handles paired-end reads and below are the analysis sections in the Snakefile:

#### &nbsp;&nbsp;**Step 1 - output files**
  - rule all - gather all output files

#### &nbsp;&nbsp;**Step 2 - gather genome data**
  - gather_genome_data: aggregate genome data from the snpeff folder
  - gatk_genome_dict: create genome dictionary
  - samtools_index: index the genome fasta file
  - bedops_gff2bed: convert the genome annotation .gff to .bed file

#### &nbsp;&nbsp;**Step 3 - quality control**
  -  trim_reads: trim adapters and low quality bases using trimmomatic or fastp

#### &nbsp;&nbsp;**Step 4 - map reads to genome**
 - bwa_index: generate bwa genome-index files for mapping
 - bwa_mem: map reads to genome, fixmate, convert .sam to .bam and remove artifacts
 - mark_duplicates: mark duplicate reads using gatk MarkDuplicatesSpark or samblaster

#### &nbsp;&nbsp;**Step 5 - mapping quality statistics**
  - samtools_idxstats: calculate alignment statistics based on the reference sequence
  - samtools_flagstats: calculates and summarizes various alignment statistics
  - samtools_depth: calculate the depth of coverage for each position in the genome
  - gatk_insert_size_metrics: collect insert size metrics

#### &nbsp;&nbsp;**Step 6 - variant calling**
  - gatk_haplotypecaller: call snps and indels via local re-assembly of haplotypes
  - generate_sample_name_map: generate a map of sample names to vcf files
  - gatk_genomics_db_import: merge gVCFs into one genomic database
  - gatk_genotype_gvcfs: perform joint genotyping

#### &nbsp;&nbsp;**Step 7 - variant filtering**
  - gatk_split_variants: separate snps and indels into separate vcf files
  - gatk_filter_hard: apply hard filters to snps and indels
  - gatk_merge_vcfs: merge snps and indels into one vcf file
  - gatk_filter_pass: filter out variants that do not pass the hard filters

#### &nbsp;&nbsp;**Step 8 - variant annotation and extraction**
  - snpeff_annotate_variants: variant annotation and functional effect prediction
  - gatk_variants_to_table: extract variant information into a table

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
  - `git clone https://github.com/kevin-wamae/gatk-variant-voyage.git`
- Type the following command in your terminal to navigate into the cloned directory using the command below. This will be the root directory of the project:
  - `cd gatk-variant-voyage`
  
- **_Note: All subsequent commands should be run from the root directory of this project. However, users can modify the scripts to their liking_**
 
 ---

## **Directory structure**
- Below is the default directory structure:
    - **config/**   - contains the Snakemake-configuration files
    - **input/** - input files
      - **bed/** - contains the bed file for specifying the intervals of interest
      - **fastq/** - contains the FastQ files
    - **output/** - contains numbered directorie of the output from the analysis
    - **workflow/** - contains the Snakefile and:
      - **envs/** - contains the Conda environment-configuration files
      - **scripts/** - contains the scripts used in the pipeline
```
.
├── LICENSE
├── README.md
├── config
├── input
│   ├── bed
│   └── fastq
├── output
└── workflow
    ├── Snakefile
    ├── envs
    └── scripts
```


- The pipeline uses _`global_wildcards()`_ to match FastQ sample names and mate files `input/fastq/` directory, using the naming convention below:
  - `reads_R1.fastq.gz` = first mate
  - `reads_R2.fastq.gz` = second mate
  - If you have a different naming convention (eg. [_this_](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm)), you can rename the FastQ files by executing the python script in the _workflow/scripts/_ directory:
    - `python workflow/scripts/fastq_rename.py`
  - Therefore, the user can deposit their FastQ files in the `input/fastq/` directory or edit the `config/config.yaml` file to point to the correct directory and the pipeline will automatically match the sample names and mates files

- The configuration file (`config/config.yaml`) specifies additional resources and can be modified to suit one's needs, such as:
  - Input files
  - Output directories,
  - The option to choose between tools, e.g.:
    - `fastp` or `trimmomatic` for read trimming
    - `gatk MarkDuplicatesSpark` or `samblaster` for marking duplicates
  - Other parameters, such as the number of threads to use

---

## **Running the analysis**
After navigating into the root directory of the project, run the analysis by executing the following commands in your terminal:

1. Create a conda analysis environment by running the command below in your terminal. This will create a conda environment named `variant-calling-gatk` and install [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [SnpEff](https://pcingola.github.io/SnpEff/se_introduction/) in it:
     - `conda env create --file workflow/envs/environment.yaml`
  
2. Activate the conda environment by running the command below in your terminal. **_Note:_** This needs to be done every time you exit and restart your terminal and want re-run this pipeline
     - `conda activate variant-calling-gatk`

3. Execute the shell script below to create the SnpEff. This will download the _P. falciparum_ genome data from [PlasmoDB](https://plasmodb.org/) and create a database in the **output/** directory. **_Note:_** This is an important step because the genome-FASTA and GFF files are required for read-mapping and variant calling.
     - `bash workflow/scripts/create_snpeff_db.sh`


4. Finally, execute the whole Snakemake pipeline by running the following command in your terminal, replace **2** with the number of CPUs you wish to use and also remember to change the number of threads in the `config/config.yaml` file:
     - `snakemake --use-conda --cores 2`
     - This will run the whole pipeline using 2 cores

5. If you want to run multiple jobs in parallel, you can use the `--jobs` flag to specify the number of jobs to run in parallel. For example, to run 2 jobs in parallel, run the following command:
    - `snakemake --use-conda --cores 2 --jobs 2`
    - This will run 2 jobs in parallel, each with 2 cores

6. Look through **output/** directory to view the results of the analysis

7. After the analysis is complete, you can deactivate the conda environment by running the following command to exit this conda environment:
     - `conda deactivate variant-calling-gatk`

---

## **Feedback and Issues**

Report any issues or bugs by openning an issue [here](https://github.com/kevin-wamae/gatk-variant-voyage/issues) or contact me via email (wamaekevin[at]gmail.com)
  
 
