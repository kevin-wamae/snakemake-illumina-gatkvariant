# **Snakemake workflow: variant calling using GATK best practices**

### **Table of contents**
- [**Snakemake workflow: variant calling using GATK best practices**](#snakemake-workflow-variant-calling-using-gatk-best-practices)
    - [**Table of contents**](#table-of-contents)
    - [**Motivation**](#motivation)
    - [**Project dependencies:**](#project-dependencies)
      - [     **Package management**](#-package-management)
      - [     **Workflow management**](#-workflow-management)
    - [**Where to start**](#where-to-start)
    - [**Directory structure**](#directory-structure)
    - [**Running the analysis**](#running-the-analysis)


### **Motivation**


- This repository contains a pipeline built with [Snakemake](https://snakemake.readthedocs.io/en/stable/) for variant calling using Illumina-generated sequences and is based on the [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery) for variant calling.
- Additionally, this pipeline aims to reproduce a recently published pipeline that optimized the GATK4 variant calling pipeline for _Plasmodium falciparum_ ([_preprint_](10.21203/rs.3.rs-2561857/v1)). However, this is not limited to _P. falciparum_ and can be used for any organism of interest.
- _**Note: The pipeline implements VCF hard-filtering, instead of the recommended soft-filtering via Variant Quality Score Recalibration (VQSR), which will be implemented in a future release**_.


- The pipeline handles paired-end reads and uses the following tools:
  - _`fastp`_ - for trimming and filtering reads
  - _`bwa`_ - for read-mapping
  - _`gatk`_ - for 'cleaning' BAM files, sorting, marking duplicates, and generating gVCFs and for variant calling and finally for filtering variants
  - `samtools` and `gatk` - for generating mapping statistics

  
- The configuration file (`config/config.yaml`) specifies additional resources and can be modified to suit one's needs, such as:
  - Input files
  - Output directories, and
  - Other parameters, such as the number of threads to use

- The pipeline also uses _`global_wildcards()`_ to match sample names and mates files in FastQ files present in the `input/`:
  - `reads_R1.fastq.gz` = first mate
  - `reads_R2.fastq.gz` = second mate
  - If the user has a different naming convention (eg. [_this_](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm)), they can rename the files using the `fastq_rename.py` script in the `scripts/` directory:
    - `python workflow/scripts/fastq_rename.py`
  - Therefore, the user can deposit their FastQ files in the `input/fastq/` directory or edit the `config/config.yaml` file to point to the correct directory and the pipeline will automatically match the sample names and mates files

---

### **Project dependencies:**

#### &nbsp;&nbsp;&nbsp;&nbsp; **Package management**
- [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) - an open-source package management system and environment management system that runs on various platforms, including Windows, MacOS, Linux


#### &nbsp;&nbsp;&nbsp;&nbsp; **Workflow management**
- [snakemake](https://github.com/snakemake/snakemake) - a workflow management system that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment, together with a clean and modern specification language in python style.

---

### **Where to start**

- Install conda for your operating System:
  - [Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
  - [MacOS](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html)
- Clone this project using the following command in your terminal:
  - `git clone https://github.com/kevin-wamae/gatk-variant-calling-hard-filter.git`
- Type the following command in your terminal to navigate into the cloned directory using the command below. This will be the root directory of the project:
  - `cd gatk-variant-calling-hard-filter`
  
- **_Note: All subsequent commands should be run from the root directory of this project. However, users can modify the scripts to their liking_**
 
 ---

### **Directory structure**
- Below is the default directory structure:
    - **config/**   - contains the Conda environment-configuration files
    - **input/** - input files
      - **bed/** - contains the bed file specifying the regions of interest
      - **fastq/** - contains the FastQ files
    - **output/** - contains numbered directorie of the output from the analysis
    - **env/**   - contains the Conda environment-configuration file
    - **scripts/** - contains the scripts used in the pipeline
    - **workflow/** - contains the Snakefile and:
      - **env/** - contains the Conda environment-configuration files
      - **scripts/** - contains the scripts used in the pipeline
```
.
├── LICENSE
├── README.md
├── config
│   └── config.yaml
├── input
│   ├── bed
│   │   └── Pf3D7_core_genome.bed
│   └── fastq
│       ├── reads_R1.fastq.gz
│       └── reads_R2.fastq.gz
├── output
└── workflow
    ├── Snakefile
    ├── env
    │   ├── bedops.yaml
    │   ├── bwa.yaml
    │   ├── environment.yml
    │   ├── gatk.yaml
    │   ├── samtools.yaml
    │   └── trimmomatic.yaml
    └── scripts
        ├── create_snpeff_db.sh
        ├── fastq_rename.py
        └── generate_sample_vcf_map.py
```

---

### **Running the analysis**
After navigating into the root directory of the project, run the analysis by executing the following commands in your terminal:

1 - Create a conda analysis environment by running the command below in your terminal. This will create a conda environment named `variant-calling-gatk` and install [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [SnpEff](https://pcingola.github.io/SnpEff/se_introduction/) in it:
  - `conda env create --file workflow/env/variant-calling-gatk.yml`
  
2 - Activate the conda environment by running the command below in your terminal. **_Note:_** This needs to be done every time you exit and restart your terminal and want re-run this pipeline
  - `conda activate variant-calling-gatk`

3 - Execute the shell script below to create the SnpEff and SnpSift database. This will download the _P. falciparum_ genome data from [PlasmoDB](https://plasmodb.org/) and create a database in the **output/** directory. **_Note:_** This is an important step because the genome-FASTA and GFF files are required for read-mapping and variant calling.
  - `bash workflow/scripts/gather_genome_files.sh`


3 - Finally, execute the whole Snakemake pipeline by running the following command in your terminal, replace **4** with the number of CPUs you wish to use and also remember to change the number of threads in the `config/config.yaml` file:
  - `snakemake --use-conda --cores 4`

4 - Look through the output files in the **output/** directory to see the results of the analysis

5 - After the analysis is complete, you can deactivate the conda environment by running the following command to exit this conda environment:
  - `conda deactivate variant-calling-gatk`


**Report any issues or bugs by openning an issue [here](https://github.com/kevin-wamae/gatk-variant-calling-for-amplicons/issues) or contact me via email (wamaekevin[at]gmail.com)**
  
 
