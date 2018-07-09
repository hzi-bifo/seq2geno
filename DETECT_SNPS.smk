#####
### Use snakemake wrappers and other mdulization method to allow the users to apply different softwares or versions
### Problems to solve: 
### 1. different softwares for the same purpose may need different options and files, such as the index files. How to specifiy them in a rule?   
#
rule process_vcf:
    input:
        vcf="{TMP_D}/{strains}/{mapper}.vcf"
    output:
        vcf_gz="{TMP_D}/{strains}/{mapper}.vcf.gz"
    shell:
        "bgzip {input[vcf]}; "
        "tabix -p vcf {output[vcf_gz]}"

rule create_vcf:
    input:
        REF=REF_FA,
        REF_FA_INDEX=REF_FA+".fai",
        BAM="{TMP_D}/{strain}/{mapper}.sorted.bam",
        BAM_INDEX="{TMP_D}/{strain}/{mapper}.sorted.bam.bai"
    output:
        vcf="{TMP_D}/{strain}/{mapper}.vcf"
    params: 
        CORES=CORES
    shell:
        "freebayes-parallel <(fasta_generate_regions.py {input[REF_FA_INDEX]} 100000) {params[CORES]} -f {input[REF]} {input[BAM]} > {output[vcf]};"

rule sort_bam:
    input:
        "{TMP_D}/{strain}/{mapper}.raw.bam"
    output:
        temp("{TMP_D}/{strain}/{mapper}.sorted.bam"),
        temp("{TMP_D}/{strain}/{mapper}.sorted.bam.bai")
    shell:
        'bamtools sort -in {input} -out {output[0]};'
        "samtools index {output[0]};"

rule mapping:
    input:
        REF=REF_FA,
        READS1=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads1'],
        READS2=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads2'],
        REF_INDEX=REF_FA+".bwt"
    output:
        temp("{TMP_D}/{strain}/{mapper}.raw.bam")
    shell:
        "bwa mem -v 2 -M -t {CORES} {input[REF]} {input[READS1]} {input[READS2]}| samtools view -b -@ {CORES} > {output[0]} ;"

rule load_reference:
    input:
        REF_FA=REF_FA
    output:
        temp(REF_FA+".fai"),
        temp(REF_FA+".bwt")
    shell:
        "samtools faidx {input[REF_FA]};" 
        "bwa index {input}"
