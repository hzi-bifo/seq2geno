# Seq2Geno

### What is Seq2Geno?
Seq2Geno is a pipeline tool for genomic analysis about bacterial populations. Compared to shell scripts, the tool test the outputs using Snakemake and manage the computational environments using Conda, which improves the reproducibility of the computational results. The output can be used as inputs of Geno2Pheno, which trains and tests phenotype predictors based on the genomic features. We aim to help the users to focus on interpreting the NGS results and gaining the biological insights.

### What does Seq2Geno cover?
- detect single nucleotide variants
- create de novo assemblies
- compute gene presence/absence and the indels
- count gene expression levels
- infer population phylogeny
- find differentially expressed genes
- conduct ancestral reconstruction of expression levels

### Get starting
- Installation of prerequisites

    - conda (tested version: 4.5.11)
    - python

- Installation of Seq2Geno

    - Download Seq2Geno:
	`git clone --recursive [project url]`
	Please don't forget the flag `--recursive`. It will help to download the required submodules (i.e. the external repositories).

    - Install the environments:
        We have packaged the required software and modules. More specifically, they are described in a yaml files that will allow conda to find suitable software and versions and install them locally (ie. requiring no `sudo`). To install them, please do:
	`conda env -f [env_file]`

    - Make it easier to use:
        We recommend to have the software searchable by editing the `.profile`:
        `echo 'export PATH=[/where/Seq2Geno/is/located/]:$PATH' >> .profile`

- Usage
Please check [wiki]

### License
Please read [the license file]

### Contact
Please open an issue if the problem cannot be solved. It will be helpful to tell us how to reproduce your problem.
