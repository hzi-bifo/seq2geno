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
# the values are passed to here and should not be seen directly by the 
# subsequent processes
user_opt= config
config= None

#####
# samples
# the variable 'STRAINS' is deprecated because the dna and rna data may include
# different samples
sys.path.insert(0, os.path.join(user_opt['seq2geno_lib_dir'],
'ParseSamplesTab'))
import ParseSamplesTab as pst
# dna
DNA_READS= pst.read_sampletab(user_opt['dna_reads'])
# rna
RNA_READS= pst.read_sampletab(user_opt['rna_reads'])

#####
# reference
REF_FA=user_opt['ref_fa']
REF_GBK=user_opt['ref_gbk']

#####
# the other paths
LIB_D= user_opt['seq2geno_lib_dir']
TMP_D='seq2geno_temp'
CORES=user_opt['cores']
STAMPY_EXE=(os.path.join(user_opt['seq2geno_lib_dir'], 
'stampy-1.0.23','stampy.py') if user_opt['stampy_exe'] is None else
user_opt['stampy_exe'])
RAXML_EXE=('raxmlHPC-PTHREADS-SSE3' if user_opt['raxml_exe'] is None else
user_opt['raxml_exe'])
RESULT_D='seq2geno'
SOFTWARE= {}
SOFTWARE['annotator']= 'prokka'
SOFTWARE['gene_sorter']= 'roary'
SOFTWARE['epr_quantifior']= 'salmon'

#####
# recruit recipes
RULE_LIB_DIR='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged/smk/rules'
RECIPE_LIB_DIR='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged/smk/recipes'

recipes= ["ng_INFER_TREE.smk", "ng_MAKE_CONS.smk", 
    "ng_DETECT_VARS.smk",
    "ng_PROCESS_VCF.smk", "ng_MASK_VCF.smk",
    "ng_CREATE_SNPS_TABLE.smk", "ng_COMPRESS_FEAT_TABLE.smk", 
    "ng_CREATE_EXPR_TABLE.smk", "LOAD_REFERENCE.smk"]

for r in recipes:
    include: os.path.join(RECIPE_LIB_DIR, r)

##### 
# determine outputs
rule all:
    input:
        os.path.join(TMP_D, 'fastAnc') if user_opt['c_ancrec'] else ''
        user_opt['tree'],
        user_opt['expr_table'],
        user_opt['syn_snps_table'],
        user_opt['syn_snps_table']+'_GROUPS',
        user_opt['syn_snps_table']+'_NON-RDNT',
        user_opt['nonsyn_snps_table'],
        user_opt['nonsyn_snps_table']+'_GROUPS',
        user_opt['nonsyn_snps_table']+'_NON-RDNT',
