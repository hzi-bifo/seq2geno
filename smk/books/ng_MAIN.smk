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
import random
import string

def random_filename():
    '''
    For generating an unique filename that fills in the rules, although it won't
be really created. 
    '''
    n= 12
    f= ''.join(random.choice(string.ascii_uppercase + string.digits) for x in
range(n))
    return(f)

#####
# config
# the config file is already one of the api argument
# No subsequent rule is allowed to access the config variables
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
# phenotypes
PHE_TABLE_F= user_opt['phe_table']

#####
# reference
REF_FA=user_opt['ref_fa']
REF_GBK=user_opt['ref_gbk']

#####
# The paths below are determined by main instead of the user
RULE_LIB_DIR=os.path.join(user_opt['seq2geno_smk_dir'], 'rules')
RECIPE_LIB_DIR=os.path.join(user_opt['seq2geno_smk_dir'], 'recipes')

#####
# command parameters
CORES=user_opt['cores']
DIFXPR_ALPHA=user_opt['dif_alpha']
DIFXPR_LFC=user_opt['dif_lfc']

#####
# the other paths
LIB_D= user_opt['seq2geno_lib_dir']
TMP_D='seq2geno_temp'
STAMPY_EXE=(os.path.join(user_opt['seq2geno_lib_dir'],
'stampy-1.0.23','stampy.py') if user_opt['stampy_exe'] is None else
user_opt['stampy_exe'])
SOFTWARE= {}
SOFTWARE['annotator']= 'prokka'
SOFTWARE['assembler']= 'spades'
SOFTWARE['gene_sorter']= 'roary'
SOFTWARE['epr_quantifior']= 'salmon'

#####
# add the output directory to the results
output_keys= ['tree', 'gpa_table', 'nonsyn_snps_table', 'syn_snps_table',
'expr_table', 'indel_table', 'dif_out', 'c_ac_out']
for k in output_keys:
    if not (user_opt[k] is None):
        user_opt[k]= os.path.join(user_opt['output_dir'], user_opt[k])

#####
# No rule is allowed to have null input or output
# Set the variables before the rules are included
TREE_OUT=random_filename() if user_opt['tree'] is None else user_opt['tree'] 
GPA_OUT=random_filename() if user_opt['gpa_table'] is None else user_opt['gpa_table']
NONSYN_SNPS_OUT=random_filename() if user_opt['nonsyn_snps_table'] is None else  user_opt['nonsyn_snps_table']
SYN_SNPS_OUT=random_filename() if user_opt['syn_snps_table'] is None else  user_opt['syn_snps_table']
EXPR_OUT=random_filename() if user_opt['expr_table'] is None else  user_opt['expr_table']
INDEL_OUT=random_filename() if user_opt['indel_table'] is None else  user_opt['indel_table']
DIF_XPR_OUT=random_filename() if user_opt['dif_out'] is None else  user_opt['dif_out']
C_ANCREC_OUT=random_filename() if user_opt['c_ac_out'] is None else  user_opt['c_ac_out']

#####
# loading all the recipes
recipes= ["LOAD_REFERENCE.smk","CONSTRUCT_ASSEMBLY.smk",
    "CREATE_INDEL_TABLE.smk",
    "CREATE_GPA_TABLE.smk", "COUNT_GPA.smk", 
    "ng_INFER_TREE.smk", "ng_MAKE_CONS.smk",
    "ng_DETECT_VARS.smk","ng_PROCESS_VCF.smk",
    "ng_MASK_VCF.smk","ng_CREATE_EXPR_TABLE.smk",
    "ng_CREATE_SNPS_TABLE.smk", "ng_COMPRESS_FEAT_TABLE.smk",
    "DIF_XPR_ANALYSIS.smk","CONT_ANC_RECONS.smk"]

for r in recipes:
    include: os.path.join(RECIPE_LIB_DIR, r)

#####
# Determine the outputs to compute
possible_targets= [user_opt['tree'], user_opt['gpa_table'],
    user_opt['nonsyn_snps_table'], user_opt['syn_snps_table'],
    user_opt['expr_table'],user_opt['indel_table']]
targets= [ f for f in possible_targets if not(f is None)]

if user_opt['dif']:
    targets.append(user_opt['dif_out'])
if user_opt['c_ac']:
    targets.append(user_opt['c_ac_out'])
if user_opt['cmpr']:
    binary_outputs= [user_opt['gpa_table'],user_opt['indel_table'],
user_opt['nonsyn_snps_table'], user_opt['syn_snps_table']]
    targets= targets+[f+'_NON-RDNT' for f in binary_outputs if not(f is None)]

#####
# lauch the workflow
rule all:
    input: targets
