#!/bin/bash
#$ -V
#$ -l h_core=20
#$ -l h_vmem=100G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v7
source activate ng_seq2geno
./bin/main.py --keep_temp --cores 20 \
  --use-ng --ref-fa data/reference/RefCln_UCBPP-PA14.fa \
  --ref-gbk data/reference/RefCln_UCBPP-PA14.gbk \
  --dna-reads /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v7/data/samples.3.dna.tsv \
  --rna-reads /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v7/data/samples.3.rna.tsv
