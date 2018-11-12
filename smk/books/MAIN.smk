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

def assign_internal_variables(user_value):
    n= 12
    f= ''.join(random.choice(string.ascii_uppercase + string.digits) for x in
range(n))
    return(f if user_value is None else user_value)

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
#sys.path.insert(0, os.path.join(user_opt['seq2geno_lib_dir'],
#'ParseSamplesTab'))
#import ParseSamplesTab as pst
sys.path.insert(0, user_opt['seq2geno_lib_dir'])
import InternalVar as iv
main_vars= iv.reader_for_main(user_opt)
main_vars.load_materials()
# dna
#DNA_READS= pst.read_sampletab(user_opt['dna_reads'])
DNA_READS= main_vars.materials.call_dna_reads()
# rna
#RNA_READS= pst.read_sampletab(user_opt['rna_reads'])
RNA_READS= main_vars.materials.call_rna_reads()
# phenotypes
#PHE_TABLE_F= user_opt['phe_table']
PHE_TABLE_F= main_vars.materials.call_phenotype_file()
print(main_vars.materials)

#####
# reference
main_vars.load_reference()
#REF_FA=user_opt['ref_fa']
REF_FA=main_vars.reference.call_ref_seq_file()
#REF_GBK=user_opt['ref_gbk']
REF_GBK=main_vars.reference.call_ref_anno_file()
print(main_vars.reference)

#####
# The paths below are determined by main instead of the user
main_vars.load_library()
#RULE_LIB_DIR=os.path.join(user_opt['seq2geno_smk_dir'], 'rules')
RULE_LIB_DIR=main_vars.library.find_rules()
#RECIPE_LIB_DIR=os.path.join(user_opt['seq2geno_smk_dir'], 'recipes')
RECIPE_LIB_DIR=main_vars.library.find_recipes()
print(main_vars.library)

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

import ExternalSoftware as es
software_pool= es.SoftwarePool(env_sensitive= True)
SOFTWARE= {}
#SOFTWARE['annotator']= 'prokka'
#SOFTWARE['assembler']= 'spades'
#SOFTWARE['gene_sorter']= 'roary'
SOFTWARE['annotator']= software_pool.find_software('prokka')
SOFTWARE['assembler']= software_pool.find_software('spades')
SOFTWARE['gene_sorter']= software_pool.find_software('roary', target_dir=
os.path.join(LIB_D, 'roary-3.8.2'))
print(SOFTWARE)

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
OUT_DIR= user_opt['output_dir']
output_opts= ['tree', 'gpa_table', 'nonsyn_snps_table', 'syn_snps_table',
    'expr_table', 'indel_table', 'dif_out', 'c_ac_out']
(TREE_OUT, GPA_OUT, NONSYN_SNPS_OUT, SYN_SNPS_OUT, EXPR_OUT, INDEL_OUT,
DIF_XPR_OUT, C_ANCREC_OUT)= [assign_internal_variables(user_opt[opt]) for opt
in output_opts]

#####
# loading all the recipes
recipes= ["LOAD_REFERENCE.smk",
    "CREATE_INDEL_TABLE.smk", "ng_COMPRESS_FEAT_TABLE.smk",
    "CREATE_SNPS_TABLE.smk", "CREATE_EXPR_TABLE.smk",
    "CREATE_GPA_TABLE.smk", "COUNT_GPA.smk", 
    "CONSTRUCT_ASSEMBLY.smk", "MAKE_CONS.smk",
    "INFER_TREE.smk","DETECT_SNPS.smk",
    "DETECT_SNPS_FOR_TABLE.smk",
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
#rule all:
#    input: targets

