#!/bin/bash
#$ -V
#$ -l h_core=16
#$ -l h_vmem=80G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v7
source activate ng_seq2geno

snakemake --notemp --snakefile=test.smk
