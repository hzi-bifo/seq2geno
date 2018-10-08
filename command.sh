#!/bin/bash
#$ -V
#$ -l h_core=20
#$ -l h_vmem=100G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v7
source activate seq2geno
./bin/main.py \
--keep_temp \
--dry-run \
--use-ng \
-cac \
  --cores 20  \
  -dx \
--ref-fa data/reference/RefCln_UCBPP-PA14.fa \
--ref-gbk data/reference/RefCln_UCBPP-PA14.gbk \
--dna-reads /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v7/data/samples.3.dna.tsv \
--rna-reads /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v7/data/samples.3.rna.tsv \
--expr results/expr.edit.tab \
--phe data/strains/pheno/phenotypes.edit.mat \
--tree results/Paeru.nwk
