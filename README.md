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
    - git (tested version: 2.17)

- Installation of Seq2Geno

    1. Download Seq2Geno:

	`git clone --recurse-submodules https://github.com/hzi-bifo/seq2geno.git`

	The flag `--recurse-submodules` will help to download the required submodules (i.e. the external repositories). The flag is available only in git version >2.13, so the users of earlier versions should considering finding the proper method to download the submodules. 

    2. Install the core environment:

	The core environment is required to initiate Seq2Geno. Please follow the few describtion in install/INSTALL

    3. Install the process-spcific environment:
	
	All the processes in Seq2Geno do not share the same pool of computational tools. The process-specific tools, however, do not need to be installed manually, because they are already listed in yaml files that Conda can parse. The installation using Conda is automatically launched when Seq2Geno is used for the first time. 

### Usage and Input

```
  seq2geno \
      -dryrun -s -d -e -p -de -ar \
      --cores 15 \
      --dna-reads samples.10.dna_full.with_ref.tsv \
      --rna-reads samples.10.rna_full.tsv \
      --ref-fa ref.fasta \
      --ref-gff ref.gff3 \
      --ref-gbk ref.gb \
      --pheno phenotype_list \
      --adaptor adapter.fasta \
      --wd SEQ2GENO_OUT
```
| option | function |
| --- | --- |
| dryrun | display the processes |
| s | SNPs calling |
| d | creating _de novo_ assemblies |
| e | counting expression levels |
| p | inferring the phylogeny |
| de | differential expression |
| ar | ancestral reconstruction of expression levels |

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

- adaptor: The adaptor file

The fasta file of adaptors of DNA-seq. It is used to process the DNA-seq reads. 

### Connection with Geno2Pheno
Geno2Pheno requires an "genml" file to describe the required files and parameters. To create the file, we offer scripts to create the file. 

  1. If your data were computed by Seq2Geno:

Please use create_genml_from_seq2geno.py
```
create_genml_from_seq2geno.py \
  --seq2geno SEQ2GENO_OUT \
  --pred abr --opt scores_f1_1 \
  --fold_n 5 --test_perc 0.1 \
  --part rand --models svm lr --k-mer 6 \
  --cls geno2pheno.classes \
  --out geno2pheno
```
| option | function |
| --- | --- |
| pred | the name of machine learning project |
| opt | the target metric to optimize |
| fold_n | the number of fold for cross validation |
| test_perc | the precentage of samples separated from the whole set for individual testing |
| part | the method to partition samples |
| models | the machine learning algorithm |
| k-mer | the k-mer size for encoding genome sequences |

- cls: classification labels
The file specifies classification labels based on the phenotypes. It is a two-column tab-separated file, where the first column are the phenotypes and the second column includes their prediction labels. 
```
high	class_1
mediate	class_1
low	class_2
```

- out: the prediction output folder

  1. If your data were not computed by Seq2Geno:

Please use create_genml.py, which allows you to specify the paths to the precomputed results. 

### License
Please read [the license file]

### Contact
Please open an issue if the problem cannot be solved. 
We will need to know how to reproduce your problem.
