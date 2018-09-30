#####
### Use snakemake wrappers and other mdulization method to allow the users to apply different softwares or versions
### Problems to solve: 
### 1. different softwares for the same purpose may need different options and files, such as the index files. How to specifiy them in a rule?   
#

rule index_vcf:
    input:
        vcf_gz="{TMP_D}/{strain}/freebayes/vcf.gz"
    output:
        vcf_gz_index= "{TMP_D}/{strain}/freebayes/vcf.gz.tbi"
    shell:
        """
        tabix -p vcf {input[vcf_gz]}
        """

rule create_vcf:
    input:
        REF=REF_FA,
        REF_FA_INDEX=REF_FA+".fai",
        BAM="{TMP_D}/{strain}/bwa/tr_sorted.bam",
        BAM_INDEX="{TMP_D}/{strain}/bwa/tr_sorted.bam.bai"
    output:
        vcf_gz="{TMP_D}/{strain}/freebayes/vcf.gz"
    params: 
        CORES=CORES
    shell:
        """
        freebayes-parallel <(fasta_generate_regions.py \
        {input[REF_FA_INDEX]} 100000) {params[CORES]} \
        -p 1 -f {input[REF]} {input[BAM]} | \
        bgzip -c > {output[vcf_gz]}
        """ 

rule sort_bam:
    input:
        paired_bam="{TMP_D}/{strain}/bwa/tr_paired.bam"
    output:
        sorted_bam=temp("{TMP_D}/{strain}/bwa/tr_sorted.bam"),
        sorted_bam_index=temp("{TMP_D}/{strain}/bwa/tr_sorted.bam.bai")
    shell:
        """
        bamtools sort -in {input} -out {output[sorted_bam]}
        samtools index {output[sorted_bam]}
        """

rule paired_mapping:
    input:
        REF=REF_FA,
        READS1=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads1'],
        READS2=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads2'],
        REF_INDEX=REF_FA+".bwt"
    output:
        paired_bam=temp("{TMP_D}/{strain}/bwa/tr_paired.bam")
    shell:
        """
        bwa mem -v 2 -M -t {CORES} {input[REF]} {input[READS1]} {input[READS2]}|\
        samtools view -b -@ {CORES} > {output[0]}
        """
