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
RESULT_D='seq2geno'
SOFTWARE= {}
SOFTWARE['annotator']= 'prokka'
SOFTWARE['gene_sorter']= 'roary'

required_smk=["ng_INFER_TREE.smk", "ng_MAKE_CONS.smk", 
"ng_DETECT_VARS.smk",
"ng_DETECT_VARS.test.smk",
"ng_PROCESS_VCF.smk", "ng_MASK_VCF.smk", "ng_CREATE_SNPS_TABLE.smk",
"ng_COMPRESS_FEAT_TABLE.smk", "ng_CREATE_EXPR_TABLE.smk", "LOAD_REFERENCE.smk"]
for smk in required_smk:
    include: os.path.join('./', smk)

rule all:
    input:
        config['tree']
        config['expr_table'],
        config['syn_snps_table'],
        config['syn_snps_table']+'_GROUPS',
        config['syn_snps_table']+'_NON-RDNT',
        config['nonsyn_snps_table'],
        config['nonsyn_snps_table']+'_GROUPS',
        config['nonsyn_snps_table']+'_NON-RDNT'
