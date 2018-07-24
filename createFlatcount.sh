#!/bin/bash
#$ -V
#$ -l h_core=10
#$ -l h_vmem=20G
cd /net/sgi/metagenomics/projects/pseudo_genomics/data/K_pneumoniae/v2
source activate phypal
parallel -j 10 --joblog coverage.log --retries 2 'samtools depth -a fastq/{}/bwa.sorted.bam > fastq/{}/bwa.coverage' ::: `ls fastq| grep -v reference`
parallel -j 10 --joblog flatcount.log --retries 2 'python bam_cov2flatcount.py -in_f fastq/{}/bwa.coverage -out_f materials/{}.flatcount' ::: `ls fastq| grep -v reference`
parallel -j 10  --joblog decompressVCF.log --retries 2 ' bgzip -d -c fastq/{}/bwa.vcf.gz > materials/{}.vcf' ::: `ls fastq| grep -v reference`

