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
        BAM="{TMP_D}/{strain}/bwa/tr_bwa.bam",
        BAM_INDEX="{TMP_D}/{strain}/bwa/tr_bwa.bam.bai"
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

rule index_bam:
    input:
        bwa_bam="{TMP_D}/{strain}/bwa/tr_bwa.bam"
    output:
        bwa_bam_index="{TMP_D}/{strain}/bwa/tr_bwa.bam.bai"
    shell:
        """
        samtools index {input.bwa_bam}
        """
rule re_redirect_bwa_result:
    input:
        bwa_raw_bam=
            lambda wildcards: "{}/{}/bwa/tr_bwa.{}.bam".format(
            wildcards.TMP_D, wildcards.strain, '1' if
            (len(DNA_READS[wildcards.strain]) == 1) else '2')
    output:
        bwa_bam="{TMP_D}/{strain}/bwa/tr_bwa.bam"
    shell:
        '''
        mv {input.bwa_raw_bam} {output.bwa_bam}
        '''

rule tr_single_mapping:
    input:
        REF=REF_FA,
        FQ1=lambda wildcards: DNA_READS[wildcards.strain][0],
        REF_INDEX=REF_FA+".bwt"
    output:
        bwa_raw_bam=temp("{TMP_D}/{strain}/bwa/tr_bwa.raw.1.bam"),
        bwa_sorted_bam=temp("{TMP_D}/{strain}/bwa/tr_bwa.sorted.1.bam"),
        bwa_nonrdnt_bam=temp("{TMP_D}/{strain}/bwa/tr_bwa.1.bam")
    params:
        cores= CORES
    shell:
        """
        bwa mem -v 2 -M -t {params.cores} {input.REF} \
{input.FQ1} |\
samtools view -b -@ {params.cores} - \
> {output.bwa_raw_bam}
        bamtools sort -in {output.bwa_raw_bam} -out {output.bwa_sorted_bam}
        samtools rmdup -s {output.bwa_sorted_bam} {output.bwa_nonrdnt_bam}
        """

rule tr_paired_mapping:
    input:
        REF=REF_FA,
        FQ1=lambda wildcards: DNA_READS[wildcards.strain][0],
        FQ2=lambda wildcards: DNA_READS[wildcards.strain][1],
        REF_INDEX=REF_FA+".bwt"
    output:
        bwa_raw_bam=temp("{TMP_D}/{strain}/bwa/tr_bwa.raw.2.bam"),
        bwa_sorted_bam=temp("{TMP_D}/{strain}/bwa/tr_bwa.sorted.2.bam"),
        bwa_nonrdnt_bam=temp("{TMP_D}/{strain}/bwa/tr_bwa.2.bam")
    params:
        cores= CORES
    shell:
        """
        bwa mem -v 2 -M -t {params.cores} {input.REF} \
{input.FQ1} {input.FQ2}|\
samtools view -b -@ {params.cores} - \
> {output.bwa_raw_bam}
        bamtools sort -in {output.bwa_raw_bam} -out {output.bwa_sorted_bam}
        samtools rmdup -S {output.bwa_sorted_bam} {output.bwa_nonrdnt_bam}
        """
