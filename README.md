# Seq2Geno

### What is it?
Seq2Geno is an automatic pipeline for detecting variants, counting gene presence/absence, and computing the transcription levels with sequencing data of bacterial population. It emphsizes reproducibility and easier accessibility. We aim to help the users to focus on interpreting the NGS results and gaining the biological insights.

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
