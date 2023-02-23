# **A Snakemake pipeline for variant calling from illumina sequences**

### **Table of contents**
- [**A Snakemake pipeline for variant calling from illumina sequences**](#a-snakemake-pipeline-for-variant-calling-from-illumina-sequences)
    - [**Table of contents**](#table-of-contents)
    - [**Motivation**](#motivation)
    - [**Project dependencies:**](#project-dependencies)
      - [     **Package management**](#-package-management)
      - [     **Workflow management**](#-workflow-management)
      - [     **Bioinformatics tools (packages)**](#-bioinformatics-tools-packages)
    - [**Where to start**](#where-to-start)
    - [**Directory structure**](#directory-structure)
    - [**Running the analysis**](#running-the-analysis)
    - [**Expected output**](#expected-output)


### **Motivation**


- This repository contains a pipeline built with Snakemake for variant calling using Illumina-generated sequences.
- It's designed for short-read sequences from _Plasmodium falciparum_ (specifically, _ama1_ and _k13_ genes). However, it can be modified to work for any other genomic loci or organism(s).


- The pipeline handles paired-end reads and uses the following tools:
  - _`fastp`_ for trimming and filtering reads
  - _`bwa`_ for read mapping
  - _`samtools`_ for BAM processing
  - _`bcftools`_ for variant calling
  - _`snpeff`_ - for annotating variants in VCFs
  - _`snpsift`_ - for extracting variants from VCFs

  
- The configuration file (`workflow/config.yaml`) specifies additional resources, parameters, input and input, but it can be modified to suit one's needs, such as:
  - input files
  - output directories, and
  - parameters such as the number of threads to use

- The pipeline also uses `global_wildcards()` to match sample names and mates files in FastQ files:
  - R1=first mates
  - R2=second mates

- This pipeline assumes that you have checked the quality of your reads using a tool such as [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). If you haven't, you can start with this earlier pipeline: [illuminaTrimWizard](https://github.com/kevin-wamae/illuminaTrimWizard)

---

### **Project dependencies:**

#### &nbsp;&nbsp;&nbsp;&nbsp; **Package management**
- [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) - an open-source package management system and environment management system that runs on various platforms, including Windows, MacOS, Linux
- _Note_: Here we are going to use `conda` but, there are other alternatives such as [_`brew`_](https://brew.sh/) 

#### &nbsp;&nbsp;&nbsp;&nbsp; **Workflow management**
- [snakemake](https://github.com/snakemake/snakemake) - a workflow management system that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment, together with a clean and modern specification language in python style.

#### &nbsp;&nbsp;&nbsp;&nbsp; **Bioinformatics tools (packages)**
- [fastp](https://github.com/OpenGene/fastp) - a tool designed to provide fast all-in-one preprocessing for FastQ files.
- [bwa](https://github.com/lh3/bwa) - an aligner for short-reads (see [_minimap2_](https://github.com/lh3/minimap2) for long-read alignment)
- [samtools](https://github.com/samtools/samtools) - a suite of programs for interacting with high-throughput sequencing data
- [bcftools](https://samtools.github.io/bcftools/bcftools.html) - a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF.
- [snpEff](https://pcingola.github.io/SnpEff/) - a genetic variant annotation and functional effect prediction toolbox.
- [snpsift](https://pcingola.github.io/SnpEff/) - a toolbox for annotating genomic variants using databases, filters, and manipulates genomic annotated variants. Once you've annotated your files using SnpEff, you can use SnpSift to help you filter large genomic datasets.

---

### **Where to start**

- Install conda for your operating System: [Windows](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html), [MacOS](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html), [Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
- Clone this project using the following command:
  - `git clone https://github.com/kevin-wamae/illuminaVariantQuest.git`
- Navigate into the cloned directory using the command below. This will be the root directory of the project:
  - `cd illuminaVariantQuest`
- **_Note: All subsequent commands should be run from the root directory of this project. However, users can modify the scripts to their liking_**
 
 ---

### **Directory structure**
- Below is the default directory structure:
    - **env/**   - contains the Conda environment-configuration file
    - **input/** - contains FastQ files
      - _Note_: There are two FastQ files here for testing, but you can add your own
    - **scripts/** - contains the scripts used in the pipeline
    - **output/** - contains the output from the analysis
    - **workflow/** - contains the Snakemake (`.smk`) script and the workflow configuration file (`.yaml`)

```
.
|-- LICENSE
|-- README.md
|-- input
|   `-- fastq
|       |-- reads_R1.fastq.gz
|       `-- reads_R2.fastq.gz
|-- output
|-- scripts
|   |-- create_snpeff_db.sh
|   `-- fastq_rename.py
`-- workflow
    |-- config.yaml
    `-- illuminaSeqAssembly.smk
```

---

### **Running the analysis**
After navigating into the root directory of the project, you can run the analysis by executing the following commands in your terminal:

1 - Create the conda analysis environment and install the dependencies from the **env/environment.yml** file by running the following command in your terminal:
  - `conda env create --file env/environment.yml`
  
2 - Activate the conda environment:
  - `conda activate illumina-variant-quest`
  - _**Note**: This needs to be done every time you want to execute this pipeline_:

3 - Execute the following shell script to create the SnpEff and SnpSift database:
  - `bash scripts/create_snpeff_db.sh`
  - This will download _P. falciparum_ data (genome FASTA, GENE and PROTEIN sequence data and the genome annotation GFF data) from [PlasmoDB](https://plasmodb.org/) and create a database in the **output/** directory
  - This script downloads release-51 of the genome data from PlasmoDB, but can be modified to download other [releases](https://plasmodb.org/common/downloads/) 
  - This is an important step because the genome-FASTA and GFF files are also needed when running the pipeline

3 - Finally, execute the whole Snakemake pipeline by running the following command in your terminal, replace **4** with the number of CPUs you wish to use:
  - `snakemake -s workflow/illuminaSeqAssembly.smk -c4`

4 - After the analysis is complete, you can deactivate the conda environment by running the following command:
  - `conda deactivate`
  
  ---
  
  ### **Expected output**
-  Below is the expected directory structure of the **output/** directory:
   -  **1_snpeff_database/** - contains the SnpEff and SnpSift database
   -  **2_genome/** - contains the genome sequences, annotation files, index files and bed files
   -  **3_trimmed_fastq/** - contains the trimmed FastQ files
   -  **4_aligned_bam/** - contains the aligned BAM files
   -  **5_map_qual_stats/** - contains the mapping quality statistics
   -  **6_vcf_files/** - contains the VCF files
   -  **7_vcfs_annotated/** - contains the annotated VCF files
   -  **8_extracted_variants/** - contains the extracted variants

```
output/
├── 1_snpeff_database
│   ├── P.falciparum
│   │   ├── cds.fa
│   │   ├── genes.gff
│   │   ├── protein.fa
│   │   ├── sequence.Pf3D7_05_v3.bin
│   │   ├── sequence.Pf3D7_06_v3.bin
│   │   ├── sequence.Pf3D7_07_v3.bin
│   │   ├── sequence.Pf3D7_08_v3.bin
│   │   ├── sequence.Pf3D7_09_v3.bin
│   │   ├── sequence.Pf3D7_10_v3.bin
│   │   ├── sequence.Pf3D7_11_v3.bin
│   │   ├── sequence.Pf3D7_12_v3.bin
│   │   ├── sequence.Pf3D7_13_v3.bin
│   │   ├── sequence.Pf3D7_14_v3.bin
│   │   ├── sequence.bin
│   │   └── snpEffectPredictor.bin
│   └── genomes
│       └── P.falciparum.fa
├── 2_genome
│   ├── annotations.bed
│   ├── annotations.gff
│   ├── genome.fa
│   ├── genome.fa.fai
│   └── regions.bed
├── 3_trimmed_fastq
│   ├── log
│   │   ├── reads.html
│   │   └── reads.json
│   ├── reads_R1.fastq.gz
│   └── reads_R2.fastq.gz
├── 4_aligned_bam
│   ├── genomeIndex
│   │   ├── genome
│   │   ├── genome.amb
│   │   ├── genome.ann
│   │   ├── genome.bwt
│   │   ├── genome.pac
│   │   └── genome.sa
│   ├── reads.bam
│   └── reads.bam.bai
├── 5_map_qual_stats
│   ├── reads.bam.flagstats.txt
│   └── reads.bam.idxstats.txt
├── 6_vcf_files
│   └── reads.vcf.gz
├── 7_vcfs_annotated
│   └── reads.vcf.gz
└── 8_extracted_variants
    ├── reads.allele.freq.txt
    └── reads.allele.txt
```
