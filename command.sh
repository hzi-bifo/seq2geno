#!/bin/bash
#$ -V
#$ -l h_core=16
#$ -l h_vmem=60G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v5
source activate seq2geno

snakemake --unlock --snakefile=MAIN.smk
snakemake --notemp --snakefile=MAIN.smk
