# Seq2Geno

### What is Seq2Geno?
Seq2Geno is a computational workflow for genomic analysis of bacterial populations. Seq2Geno uses Snakemake to schedule the processes and relies on Conda to manage the computational environments. The two tools improve the reproducibility of the computational results. 
Seq2Geno outputs are formatted for use with the Geno2Pheno workflow, which trains phenotype predictors based on genomic features. 

### What does Seq2Geno cover?
- detect single nucleotide variants
- create de novo assemblies
- compute gene presence/absence and the indels
- count gene expression levels
- infer population phylogeny
- find differentially expressed genes
- reconstruct ancestral values of expression level

### Get started
- Prerequisites

    - conda (tested version: 4.6.14)
    - python (tested verson: 3.6)
    - Linux (tested version: Debian GNU/Linux 8.8 jessie)
    - git (tested version: 2.18)

- Installation of Seq2Geno

    1. Download Seq2Geno:

	```
	git clone --recurse-submodules https://github.com/hzi-bifo/seq2geno.git
	cd seq2geno
	git submodule update --init --recursive
	```

	The option `--recurse-submodules` helps to download the submodules that are located at another repository (i.e. Seq2Geno and Geno2Pheno). The flag is available only in git version >2.13, and users of earlier git versions may need to find the substitute.  

    2. Install the core environment:

	The core environment is required to initiate Seq2Geno. Please follow the few describtion in install/INSTALL

    3. Install the process-spcific environment:
	
	All the processes in Seq2Geno do not share the same pool of computational tools. The process-specific tools, however, do not need to be installed manually, because they are already listed in yaml files that Conda can parse. The installation using Conda is automatically launched when Seq2Geno is used for the first time. 

### Usage and Input

Usage:
```
  seq2geno -f options.yml
```

The input file is an yaml file where all options are described. The file includes two parts:

1. features:

| option | action | values ([default])|
| --- | --- | --- |
| dryrun | display the processes and exit | [Y]/N |
| snps | SNPs calling | Y/[N] |
| denovo | creating _de novo_ assemblies | Y/[N] |
| expr | counting expression levels | Y/[N] |
| phylo | inferring the phylogeny | Y/[N] |
| de | differential expression | Y/[N] |
| ar | ancestral reconstruction of expression levels | Y/[N] |

2. general: 

    - cores: accessible number of cpus

    - old_config: do not overwrite the previous config files. Please make the old config file really are accessible.

    - wd: the working directory
    The intermediate and final files will be stored under this folder. The final outcomes will be symlinked to RESULTS/.

    - dna_reads: The list of DNA-seq data 

    It should be a two-column list, where the first column includes all samples and the second column lists the __paired-end reads files__. The two reads file are separated by a comma. The first line is the first sample.
    ```
    sample01	/paired/end/reads/sample01_1.fastq.gz,/paired/end/reads/sample01_2.fastq.gz
    sample02	/paired/end/reads/sample02_1.fastq.gz,/paired/end/reads/sample02_2.fastq.gz
    sample03	/paired/end/reads/sample03_1.fastq.gz,/paired/end/reads/sample03_2.fastq.gz
    ```

    - rna_reads: The list of RNA-seq data

    It should be a two-column list, where the first column includes all samples and the second column lists the __short reads files__. The first line is the first sample.
    ```
    sample01	/transcription/reads/sample01.rna.fastq.gz
    sample02	/transcription/reads/sample02.rna.fastq.gz
    sample03	/transcription/reads/sample03.rna.fastq.gz
    ```

    - phe_table: The phenotype table

    The table is tab-separated. For n samples with m phenotypes, the table is (n+1)-by-(m+1) as shown below. The first column should be sample names. The header line should includes names of phenotypes. Missing values are acceptable.
    ```
    strains	virulence
    sample01	high
    sample02	mediate
    sample03	low
    ```

    - ref_fa, ref_gff, ref_gbk: the data of reference genome

    The fasta, gff, and genbank files of a reference genome. They should have same sequence ids. 

    - adaptor: The adaptor file (optional)

    The fasta file of adaptors of DNA-seq. It is used to process the DNA-seq reads. 

### License
Apache 2.0 (please see the LICENSE)

### Contact
Please open an issue here or contact Tzu-Hao Kuo (Tzu-Hao.Kuo@helmhotz-hzi.de). 
We will need to know how to reproduce the problem. 

