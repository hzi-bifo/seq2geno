import pandas as pd
import os
import re
configfile: "config.yaml"
SAMPLES_DF=pd.read_table(config["samples"], sep= '\t', header= 0).set_index("strain", drop=False)
STRAINS=SAMPLES_DF['strain'].tolist()
REF_FA=config['reference_sequence']
REF_GBK=config['reference_annotation']
TMP_D=(config['tmp_d'] if re.search('\w', config['tmp_d']) else '.')
CORES=config['cores']
STAMPY_EXE=config['stampy_exe']
RAXML_EXE=config['raxml_exe']
TREE=config['tree']
GPA_TABLE=config['gpa_table']
RESULT_D=config['result_d']

include: "CREATE_GPA_TABLE.smk"
include: "COUNT_GPA.smk"
include: "CONSTRUCT_ASSEMBLY.smk"
include: "MAKE_CONS.smk"
include: "INFER_TREE.smk"
include: "DETECT_SNPS.smk"

rule all:
    input:
        #vcf_gz=expand(TMP_D+"/{strains}/{mapper}.vcf.gz", strains= STRAINS, mapper= 'bwa')
        #cons_seqs=expand(TMP_D+"/{strains}/{mapper}.cons.fa", strains= STRAINS, mapper= 'bwa')
        TREE,
#        ASSEM=expand(TMP_D+"/{strains}/{assembler}.assem.fa", strains= STRAINS, assembler= 'spades'),
#        roary_gpa=TMP_D+"/roary/gene_presence_absence.csv"
        GPA_TABLE
        
    
'''
rule import_dna_reads:
    input:
    output:
        dna_reads

rule import_dna_params:
    input:

    output:
        dna_params

rule import_reference:
    input:

    output:
        ref

rule import_rna_reads:
    input:
    
    output:
        rna_reads

rule import_rna_params:
    input:

    output:
        rna_params

rule count_gpa:
    input:
        assembly
    output:
        roary_gpa

rule construct_assembly:
    input:
        dna_reads
        dna_params
    output:
        assembly

rule detectSNPs:
    input:
        dna_reads
        dna_params
        ref
    output:
        vcf

rule compute_expressions:
    input:
        rna_reads
        rna_params
        ref
    output:

rule create_gpa_table:
    input:
        roary_gpa
    output:
        gpa_table

rule create_snps_table:
    input:
        vcf
    output:
        snps_table

rule infer_tree:
    input:
        one_big_aln
    output:
        tree

rule make_cons_seqs:
    input:
        vcf
        ref
    output:
        cons_seqs=[...]



rule create_expr_table:
    input:

    output:
        expr_table
             
rule make_ml_input:
    input:
        gpa_table
        snps_table
        tree
        expr_table
        assembly
    output:
        ml_markdown
'''
