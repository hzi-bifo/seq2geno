#!/bin/bash
#$ -V
#$ -l h_core=20
#$ -l h_vmem=100G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged
#source activate seq2geno
./bin/main.py \
  --ng\
  --c_ac \
  --dx \
  --dryrun \
  --keeptemp \
  --cores 20  \
  --ref-fa data/reference/RefCln_UCBPP-PA14.fa \
  --ref-gbk data/reference/RefCln_UCBPP-PA14.gbk \
  --dna-reads /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged/data/samples.3.dna.tsv \
  --rna-reads /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged/data/samples.3.rna.tsv \
  --expr results/expr.edit.tab \
  --phe data/strains/pheno/phenotypes.edit.mat \
  --tree results/Paeru.nwk
