#!/bin/bash
#$ -V
#$ -l h_core=15
#$ -l h_vmem=200G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno.dirty.in_one
source activate py37
source activate snakemake_env
export SEQ2GENO_HOME=/net/metagenomics/data/from_moni/old.tzuhao

main/seq2geno \
  --dryrun \
  --cores 15 \
  --dna-reads /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/data/samples.10.dna_full.with_ref.tsv \
  --rna-reads /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/data/samples.10.rna_full.tsv \
  --ref-fa /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/data/reference/Pseudomonas_aeruginosa_PA14.edit.fasta \
  --ref-gff /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v9/.ori_backup/gpa/match_reference_annotation/RefCln_UCBPP-PA14.edit.gff \
  --ref-gbk /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/data/reference/Pseudomonas_aeruginosa_PA14_ncRNA.gbk \
  --pheno /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/data/strains/pheno/phenotypes.edit.mat \
  --wd whatever
create_genml/create_genml_from_seq2geno.py \
  --genml ./test.gml --seq2geno whatever \
  --pred abr --opt scores_f1_1 \
  --fold_n 5 --test_perc 0.1 \
  --part rand --models svm lr --k-mer 6 \
  --cls create_genml/classes --out geno2pheno
