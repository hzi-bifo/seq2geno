# Seq2Geno

### What is Seq2Geno?
Seq2Geno is a pipeline tool for genomic analysis about bacterial populations. Compared to shell scripts, the tool tests the outputs using Snakemake and manages the computational environments using Conda. The two tools improve the reproducibility of the computational results. The outputs from Seq2Geno can be used as inputs of Geno2Pheno, which trains and tests phenotype predictors based on the genomic features. We aim to help the users to focus on interpreting the NGS results and gaining the biological insights.

### What does Seq2Geno cover?
- detect single nucleotide variants
- create de novo assemblies
- compute gene presence/absence and the indels
- count gene expression levels
- infer population phylogeny
- find differentially expressed genes
- conduct ancestral reconstruction of expression levels

### Get starting
- each process are described in yaml files, which include the software and pacakges. These yaml files will be parsed by Conda and create the computatioal environments. 
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

### Connection with Geno2Pheno
Geno2Pheno requires an "genml" file to describe the required files and parameters. To create the file, we offer scripts to create the file. 
- If your data were computed by Seq2Geno:
Please use create_genml_from_seq2geno.py
- If your data were not computed by Seq2Geno:
Please use create_genml.py, which allows you to specify the paths to the precomputed results. 

### License
Please read [the license file]

### Contact
Please open an issue if the problem cannot be solved. 
We will need to know how to reproduce your problem.
