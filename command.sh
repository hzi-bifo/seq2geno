#!/bin/bash
#$ -V
#$ -l h_core=20
#$ -l h_vmem=100G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged
#source activate seq2geno
./bin/main.py \
  --c_ac \
  --dx \
  --dryrun \
  --keeptemp \
  --cores 20  \
  --ref-fa data/reference/RefCln_UCBPP-PA14.fa \
  --ref-gbk data/reference/RefCln_UCBPP-PA14.gbk \
  --dna-reads /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged/data/samples.3.dna.tsv \
  --rna-reads /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged/data/samples.3.rna.tsv \
  --phe data/strains/pheno/phenotypes.edit.mat \
  --outdir results \
  --expr expr.edit.tab \
  --tree Paeru.nwk \
  --s-snp s-syn.tab 
