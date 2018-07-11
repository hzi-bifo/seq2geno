#!/bin/bash
#$ -V
#$ -l h_core=16
#$ -l h_vmem=60G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3
source activate phypal
snakemake  --snakefile=INTERFACE.smk
