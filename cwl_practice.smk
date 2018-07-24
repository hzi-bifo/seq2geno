import pandas as pd
import os

configfile: "config.yaml"
SAMPLES_DF=pd.read_table(config["samples"], sep= '\t', header= 0).set_index("strain", drop=False)
STRAINS=SAMPLES_DF['strain'].tolist()
REF_FA=config['reference_sequence']
REF_GBK=config['reference_annotation']
TMP_D=config['tmp_d']
RESULT_D=config['result_d']
CORES=config['cores']
STAMPY_EXE=config['stampy_exe']
RAXML_EXE=config['raxml_exe']
PREFIX= 'mix'

rule all:
    input:
        expand("{tmp_dir}/{strain}/{mapper}.sorted.bam", tmp_dir= TMP_D, strain= STRAINS, mapper= 'bwa')    

rule prepare_reference:
    input:
        REF_FA
    output:
        REF_FA+".fai",
        REF_FA+".bwt"

    shell:
        "samtools faidx {input};" 
        "bwa index {input}"

rule map2reference:
    input:
        REF=REF_FA,
        READS1=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads1'],
        READS2=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads2'],
        REF_INDEX=REF_FA+".bwt"
    output:
        temp("{tmp_dir}/{strain}/{mapper}.raw.bam")
#    shell:
#        "bwa mem -v 2 -M -t {CORES} {input[REF]} {input[READS1]} {input[READS2]}| samtools view -b -@ {CORES} > {output[0]} ;"
    run:
        if == 'bwa':
            

rule samtools_sort:
    input:
        input="{tmp_dir}/{strain}/{mapper}.raw.bam"
    output:
        output_name="{tmp_dir}/{strain}/{mapper}.sorted.bam"
    params:
        threads=lambda wildcards, threads: threads,
        runtime.tmpdir=wildcords.tmp_dir,
        memory="4G"
    threads: 8
    cwl:
        "https://github.com/common-workflow-language/workflows/blob/"
        "fb406c95/tools/samtools-sort.cwl"

