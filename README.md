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

    - conda (tested version: 4.5.11)
    - python (tested verson: 3.6)
    - Linux (tested version: Debian GNU/Linux 8.8 jessie)
    - git (tested version: 2.18)

- Installation of Seq2Geno

    1. Download Seq2Geno:

	`git clone --recurse-submodules https://github.com/hzi-bifo/seq2geno.git`

	The flag `--recurse-submodules` will help to download the required submodules (i.e. the external repositories). The flag is available only in git version >2.13, so the users of earlier versions should considering finding the proper method to download the submodules. 

    2. Install the core environment:

	The core environment is required to initiate Seq2Geno. Please follow the few describtion in install/INSTALL

    3. Install the process-spcific environment:
	
	All the processes in Seq2Geno do not share the same pool of computational tools. The process-specific tools, however, do not need to be installed manually, because they are already listed in yaml files that Conda can parse. The installation using Conda is automatically launched when Seq2Geno is used for the first time. 

### Usage and Input

Usage:
```
  seq2geno -f options.yml
```

The input file is an yaml file where all options are described. The file includes two parts. 
- functions:

| option | action | values ([default])|
| --- | --- | --- |
| dryrun | display the processes and exit | [Y]/N |
| s | SNPs calling | Y/[N] |
| d | creating _de novo_ assemblies | Y/[N] |
| e | counting expression levels | Y/[N] |
| p | inferring the phylogeny | Y/[N] |
| de | differential expression | Y/[N] |
| ar | ancestral reconstruction of expression levels | Y/[N] |

- files:

    - cores: available number of cpus 
    Although the parameter is included in the "files" session, please just set a number instead of specify an file where the number is stated.

    - wd: the working directory
    The intermediate and final files will be created under the folder. The final outcomes will be symlinked to RESULTS/.

    - dna-reads: The list of DNA-seq data 

    It should be a two-column list, where the first column includes all samples and the second column lists the __paired-end reads files__. The two reads file are separated by a comma. The first line is the first sample.
    ```
    sample01	sample01_1.fastq.gz,sample01_2.fastq.gz
    sample02	sample02_1.fastq.gz,sample02_2.fastq.gz
    sample03	sample03_1.fastq.gz,sample03_2.fastq.gz
    ```

    - rna-reads: The list of RNA-seq data

    It should be a two-column list, where the first column includes all samples and the second column lists the __short reads files__. The first line is the first sample.
    ```
    sample01	sample01.rna.fastq.gz
    sample02	sample02.rna.fastq.gz
    sample03	sample03.rna.fastq.gz
    ```

    - pheno: The phenotype table

    For n samples with m phenotypes, the table is n-by-m and the upper-left is blanck. The table is tab-separated. The header line includes the name of phenotypes. The value in the table can also be blanck. 
    ```
	    virulence
    sample01	high
    sample02	mediate
    sample03	low
    ```

    - ref-fa, ref-gff, ref-gbk	The reference data

    The fasta, gff, and genbank files of a reference genome. They should have same sequence ids. 

    - adaptor: The adaptor file (optional)

    The fasta file of adaptors of DNA-seq. It is used to process the DNA-seq reads. 

### License
Please read [the license file]

### Contact
Please open an issue here, or send an email to Tzu-Hao.Kuo@helmhotz-hzi.de

