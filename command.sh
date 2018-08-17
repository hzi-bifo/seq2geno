#!/bin/bash
#$ -V
#$ -l h_core=16
#$ -l h_vmem=60G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v4 
source activate seq2geno
snakemake --notemp --snakefile=INTERFACE.smk
