#!/bin/bash
#$ -V
#$ -l h_core=20
#$ -l h_vmem=100G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v9
#source activate seq2geno
./bin/seq2geno \
  --project test \
  --ref-fa data/reference/RefCln_UCBPP-PA14.fa \
  --ref-gbk data/reference/RefCln_UCBPP-PA14.gbk \
  --dna-reads /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged/data/samples.3.dna.tsv \
  --rna-reads /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged/data/samples.3.rna.tsv \
  --phe data/strains/pheno/phenotypes.edit.mat \
  --c_ac \
  --dx \
  --dryrun \
  --keeptemp \
  --cores 20  \
  --outdir results \
  --tree Paeru.nwk \
  --gpa gpa.mat \
  --ns-snp ns-syn.mat \
  --s-snp s-syn.mat \
  --expr expr.edit.tab \
  --ind indel.mat 
