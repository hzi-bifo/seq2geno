#####
### Use snakemake wrappers and other mdulization method to allow the users to apply different softwares or versions
### Problems to solve: 
### 1. different softwares for the same purpose may need different options and files, such as the index files. How to specifiy them in a rule?   

rule all:
    input:
        vcf_gz

rule process_vcf:
    input:
        vcf
    output:
        vcf_gz

rule create_vcf:
    input:

    output:

rule sort_bam:
    input:
    
    output:

rule mapping:
    input:

    output:

rule load_reads
    input:

    output:

rule load_reference
    input:
        ref_fa
    output:
        
'''  
rule prepare_reference:
    input:
        REF_FA
    output:
        REF_FA+".fai",
        REF_FA+".bwt"

    shell:
        "samtools faidx {input};" 
        "bwa index {input}"

rule prepare_reference_Ariane:
    input:
        REF_FA
    output:
        REF_FA+".stidx",
        REF_FA+".sthash"
    params:
        STAMPY=STAMPY_EXE
    shell:
        "source activate py27;"
        "{params[STAMPY]} -G {input} {input};"
        "{params[STAMPY]} -g {input} -H {input};"
        "source deactivate"

rule map2reference_Ariane:
    input:
        REF=REF_FA,
        READS1=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads1'],
        READS2=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads2'],
        REF_INDEX=REF_FA+".stidx"
    output:
        "{tmp_dir}/{strain}/stampy.sorted.bam",
        "{tmp_dir}/{strain}/stampy.sorted.bam.bai"
    params:
        STAMPY=STAMPY_EXE,
    	CORES=CORES
    shell:
        # stampy runs with python2.6
        "source activate py27;"
        "{params[STAMPY]} --bwaoptions=\"-q10 {input[REF]}\" -g {input[REF]} -h  {input[REF]} -M  {input[READS1]} {input[READS2]} | samtools view -bS -@ {params[CORES]} |  samtools sort>  {output[0]};"
        "samtools index {output[0]};"
        "source deactivate;"

rule call_var:
    input:
        REF=REF_FA,
        REF_FA_INDEX=REF_FA+".fai",
        BAM="{tmp_dir}/{strains}/stampy.sorted.bam",
        BAM_INDEX="{tmp_dir}/{strains}/stampy.sorted.bam.bai"
    output:
        "{tmp_dir}/{strains}/"+PREFIX+".vcf"
    params: 
        CORES=CORES
    shell:
        "freebayes-parallel <(fasta_generate_regions.py {input[REF_FA_INDEX]} 100000) {params[CORES]} -f {input[REF]} {input[BAM]} > {output};"

rule compress_vcf:
    input:
        "{tmp_dir}/{strain}/{PREFIX}.vcf"
    output:
        #TMP_D+"/{strain}/bwa.vcf.gz",
        #TMP_D+"/{strain}/bwa.vcf.gz.tbi"
        "{tmp_dir}/{strain}/{PREFIX}.vcf.gz",
        "{tmp_dir}/{strain}/{PREFIX}.vcf.gz.tbi"
    shell:
        "bgzip {input}; "
        "tabix -p vcf {output[0]}"

rule process_vcf:
    input:
        TMP_D+"/{strain}/{PREFIX}.vcf.gz"
    output:
        TMP_D+"/{strain}/{PREFIX}.vcf.non-indel.gz",
'''
