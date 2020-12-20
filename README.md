# Seq2Geno

- [What is Seq2Geno?](#intro)
- [What does Seq2Geno cover?](#functions) 
- [Get started](#install) 
- [Usage and Input](#usage) 
    - [GUI](#gui)
    - [command line](#commandline)
    - [arguments](#args)
- [Train the phenotypic predictor with the Seq2Geno results](#genyml) 
- [Example usages and data](#example) 
- [FAQ](#FAQ)
- [License](#license) 
- [Contact](#contact) 


### <a name="intro"></a>What is Seq2Geno?
Seq2Geno2Pheno complements experiemental methods to facilitate the study of genotypes and phenotypes. As the first stage of Seq2Geno2Pheno, Seq2Geno provides user-friendly access to computational analyses with the next-generation sequencing data of bacterial samples. It enables the users to use either the graphical or the command line interface to edit the arguments and conduct the complex analyses. The final output data include the phylogenetic tree and the feature matrices of SNPs, gene presence/absence, and expression levels, which are formatted for the next stage of machine learning: Geno2Pheno.

Besides the easy access, expanding the input dataset could be also effortless (See section [FAQ](#FAQ) for more details). Furthermore, Seq2Geno automatically resolves the dependencies among procedures and manages the procedure-specific computational environments, aiming to avoid error-prone behaviors of researchers (such as manually repeating processes or shifting computational environments). 

The output data from Seq2Geno can be used to train phenotypic predictors using [Geno2Pheno](https://genopheno.bifo.helmholtz-hzi.de). To have your data submitted to the Geno2Pheno server, please check [the section below](#genyml). 

### <a name="functions"></a>What does Seq2Geno cover?
- detect single nucleotide variants
- create de novo assemblies
- compute gene presence/absence and the indels
- count gene expression levels
- infer the phylogenetic tree
- find differentially expressed genes (additional data that won't be used by Geno2Pheno)
- reconstruct ancestral values of expression level (additional data that won't be used by Geno2Pheno)

### <a name="install"></a>Get started
- Prerequisites

    - [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) (tested version: 4.8.4)
    - file [.condarc](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html) that includes these channels and is detectable by your conda
      - hzi-bifo
      - conda-forge/label/broken
      - bioconda
      - conda-forge
      - defaults
    - [python](https://www.python.org/downloads/) (tested verson: 3.6)
    - [Linux](https://www.cyberciti.biz/faq/find-linux-distribution-name-version-number/) (tested version: Debian GNU/Linux 8.8 jessie)
    - [git](https://git-scm.com/downloads) (tested version: 2.18)

- Installation of Seq2Geno

	Seq2Geno doen't need compilation. Therefore, you will only need to clone Seq2Geno:
	```
	git clone --recurse-submodules https://github.com/hzi-bifo/seq2geno.git
	cd seq2geno
	git submodule update --init --recursive
	```

	The option `--recurse-submodules` helps to download the submodules that are located at another repository (i.e. Seq2Geno and Geno2Pheno). The flag is available only in git version >2.13, and users of earlier git versions may need to find the substitute.  

- Installation of the environments 

Please go to `install/` and either use the `INSTALL.sh` or follow the manual `INSTAL.md`. 

### <a name="usage"></a>Usage and Input

The graphical user interface (GUI) `seq2geno_gui` and the command line `seq2geno` can be either launched using the launcher `S2G`, put under the home folder of seq2geno, or found in the main folder. When no argument is set for `S2G`, the GUI will be launched; otherwise, it passes arguments to the command line tool:

```
  S2G -d -f [options_yaml] -l [log_file]
```

Please read the subset about [command line](#commandline) for more information.
 

To use the two interfaces without the launcher, please remeber to activate the core environment:

```
conda activate snakemake_env
``` 

or 
```
source activate snakemake_env
```

. More information about usages are listed below.

- <a name="gui"></a>GUI

Use the tool `seq2geno_gui` to read, edit, or save the arguments in a yaml file. Once the arguments are ready, the analyses can be launched with this interface; for large-scale researches, however, generating the yaml file and launching the analyses with the command line method (described below) might be more convenient, as having processes running in background should be more convenient. To learn more, please read the the manual `doc/GUI_manual.pdf`.

- <a name="commandline"></a>command line

The input for seq2geno is a single yaml file describing all arguments:
```
  seq2geno -d -f [options_yaml] -l [log_file]
```

The [options\_yaml] describes all the options and input data for Seq2Geno. The [log\_file] should be a non-existing filename to store the log information; if not set, the messages will be directed to stdout and stderr.

- <a name="args"></a>arguments

The input file is an yaml file where all options are described. The file includes two parts:

1. features:

| option | action | values ([default])|
| --- | --- | --- |
| dryrun | display the processes and exit | [Y]/N |
| snps | SNPs calling | Y/[N] |
| denovo | creating de novo assemblies | Y/[N] |
| expr | counting expression levels | Y/[N] |
| phylo | inferring the phylogeny | Y/[N] |
| de | differential expression | Y/[N] |
| ar | ancestral reconstruction of expression levels | Y/[N] |

To only create the folder and config files, please turn off the last six options. 

2. general (\* mandatory): 

    - cores: number of cpus (integer; automatically adjusted if larger than the available cpu number)

    - mem_mb: memory size to use (integer in mb; automatically adjusted if larger than the free memory). __Note: some processes may crush because of insufficiently allocated  memory__

    - \*wd: the working directory. The intermediate and final files will be stored under this folder. The final outcomes will be symlinked to the sub-directory RESULTS/.

    - \*dna_reads: the list of DNA-seq data 

    It should be a two-column list, where the first column includes all samples and the second column lists the __paired-end reads files__. The two reads file are separated by a comma. The first line is the first sample.
    ```
    sample01	/paired/end/reads/sample01_1.fastq.gz,/paired/end/reads/sample01_2.fastq.gz
    sample02	/paired/end/reads/sample02_1.fastq.gz,/paired/end/reads/sample02_2.fastq.gz
    sample03	/paired/end/reads/sample03_1.fastq.gz,/paired/end/reads/sample03_2.fastq.gz
    ```

    - \*ref_fa, ref_gff, ref_gbk: the data of reference genome

    The fasta, gff, and genbank files of a reference genome. They should have same sequence ids. 


    - old_config: if recognizable, the config files that were previously stored in the working directory will be reused. ('Y': on; 'N': off)

    - rna_reads: the list of RNA-seq data. (string of filename)

    It should be a two-column list, where the first column includes all samples and the second column lists the __short reads files__. The first line is the first sample.
    ```
    sample01	/transcription/reads/sample01.rna.fastq.gz
    sample02	/transcription/reads/sample02.rna.fastq.gz
    sample03	/transcription/reads/sample03.rna.fastq.gz
    ```

    - phe_table: the phenotype table (string of filename)

    The table is tab-separated. For n samples with m phenotypes, the table is (n+1)-by-(m+1) as shown below. The first column should be sample names. The header line should includes names of phenotypes. Missing values are acceptable.
    ```
    strains	virulence
    sample01	high
    sample02	mediate
    sample03	low
    ```

    - adaptor: the adaptor file (string of filename)

    The fasta file of adaptors of DNA-seq. It is used to process the DNA-seq reads. 

### <a name="genyml"></a>Train the phenotypic predictor with the Seq2Geno results 
[Geno2Pheno](https://genopheno.bifo.helmholtz-hzi.de) requires all input data packed in a single zip file. The input file for the validator and gnerator of that zip file can be generated using submission\_tool/create\_genyml.py. 

### <a name="example"></a>Example usages and data
The tutorials and example data and commands can be found in  `example_sg_dataset`. In case the folder was not decompressed from the tar.gz file, please

```
tar zxvf ./example_sg_dataset.tar.gz 
```
### <a name="FAQ"></a>FAQ
__Will every procedure be rerun if I want to add one sample?__

No, you will need to add one more line in your reads list (i.e., the dna or the rna reads. See section [arguments](#args) for more details.) and then run the same workflow again. Seq2Geno use Snakemake to determine whether certain intermediate data need to be recomputed or not. 

__Will every procedure be rerun if I want to exclude one sample?__

No; however, besides excluding that sample from the reads list, you will need to remove all the subsequent results that were previously computed. That could be risky.

__Will every procedure be rerun if I accidentally delete some intermediate data?__

No, only the deleted one and the subsequent data will be recomputed.

__Where is the final data?__

Under your working directory, they are collected in the subfolder `RESULTS/`.

__What is the current status?__

If the log file was specified when Seq2Geno was launched, you could check the log file to determine the current status. Otherwise, the status should be directed to your STDOUT or STDERR.

__Note: you might need to ensure the memory setting in the seq2geno_input.yml__

### <a name="license"></a>License
Apache 2.0 (please see the LICENSE)

### <a name="contact"></a>Contact
Please open an issue here or contact Tzu-Hao Kuo (Tzu-Hao.Kuo@helmhotz-hzi.de). 
We will need to know how to reproduce the problem. 

