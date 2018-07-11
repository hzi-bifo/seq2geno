#!/bin/bash
#$ -V
#$ -l h_vmem=150G

cd /net/sgi/metagenomics/projects/pseudo_genomics/data/K_pneumoniae/v2/materials.v2

python2 ./mutation_table.v4.py  -f strain_list -a anno.tab -o all_SNPs.v4.tab
#python2 ./mutation_table.py  -f strain_list -a anno.tab -o all_SNPs.tab
python2 ./Snp2Amino.py -f all_SNPs.v4.tab -g reference.gbk -n all -o syn_SNPs_final.v4.tab
python2 ./Snp2Amino.py -f all_SNPs.v4.tab -g reference.gbk -n non-syn -o non-syn_SNPs_final.v4.tab
