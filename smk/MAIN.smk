'''
It is an interface, or convertor, to pass the users requests that were already
parsed by the main script to the snakemake rules. It determines all the
snakemake workflows and the global variables required by the subsequent
processes, including those described in config and samples table. It is the
backbone of all the computatinal processes.

Input: 
    arguments from the main script

Output:
    user-interested genomic data
    
'''

import pandas as pd
import os
import re
import sys

#####
# config
# the config file is already one of the api argument

#####
# samples
# the variable 'STRAINS' is deprecated because the dna and rna data may include
# different samples
sys.path.insert(0, os.path.join(config['seq2geno_lib_dir'],
'ParseSamplesTab'))
import ParseSamplesTab as pst
# dna
DNA_READS= pst.read_sampletab(config['dna_reads'])
# rna
RNA_READS= pst.read_sampletab(config['rna_reads'])

#####
# reference
REF_FA=config['ref_fa']
REF_GBK=config['ref_gbk']

#####
# the other paths
TMP_D='seq2geno_temp'
CORES=config['cores']
STAMPY_EXE=(os.path.join(config['seq2geno_lib_dir'], 
'stampy-1.0.23','stampy.py') if config['stampy_exe'] is None else
config['stampy_exe'])
RAXML_EXE=('raxmlHPC-PTHREADS-SSE3' if config['raxml_exe'] is None else
config['raxml_exe'])


include: "LOAD_REFERENCE.smk"
include: "CREATE_INDEL_TABLE.smk"
rule all:
    input:
        expand('{TMP_D}/{strain}/stampy/dna_for_tab.sam', 
        TMP_D=TMP_D, strain= DNA_READS.index.values.tolist())
'''
include: "LOAD_REFERENCE.smk"
include: "CREATE_INDEL_TABLE.smk"
include: "CREATE_SNPS_TABLE.smk"
include: "CREATE_EXPR_TABLE.smk"
include: "CREATE_GPA_TABLE.smk"
include: "COUNT_GPA.smk"
include: "CONSTRUCT_ASSEMBLY.smk"
include: "MAKE_CONS.smk"
include: "INFER_TREE.smk"
include: "DETECT_SNPS.smk"
include: "DETECT_SNPS_FOR_TABLE.smk"

rule all:
    input:
        config['expr_table'],
        config['syn_snps_table'],
        config['syn_snps_table']+'_GROUPS',
        config['syn_snps_table']+'_NON-RDNT',
        config['nonsyn_snps_table'],
        config['nonsyn_snps_table']+'_GROUPS',
        config['nonsyn_snps_table']+'_NON-RDNT',
        config['tree'],
'''
