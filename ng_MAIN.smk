import pandas as pd
import os
import re
configfile: "config.yaml"
SAMPLES_DF=pd.read_table(config["samples"], sep= '\t', header= 0).set_index("strain", drop=False)
STRAINS=SAMPLES_DF['strain'].tolist()
REF_FA=config['reference_sequence']
REF_GBK=config['reference_annotation']
#TMP_D=(config['tmp_d'] if re.search('\w', config['tmp_d']) else '.')
TMP_D='seq2geno_temp'
CORES=config['cores']
STAMPY_EXE=config['stampy_exe']
RAXML_EXE=config['raxml_exe']
#RESULT_D=config['result_d']
RESULT_D='seq2geno'
#SOFTWARE={'mapper': 'bwa'}
SOFTWARE= config['software']
SOFTWARE['annotator']= 'prokka'
SOFTWARE['gene_sorter']= 'roary'
print(SAMPLES_DF)
print(SOFTWARE)

<<<<<<< HEAD:ng_MAIN.smk
include: "ng_INFER_TREE.smk"
include: "ng_MAKE_CONS.smk"
=======
>>>>>>> 34a99652488439fc7c9af0fa9c8ceb3cfc128396:ng_MAIN.smk
include: "ng_DETECT_VARS.smk"
include: "ng_PROCESS_VCF.smk"
include: "ng_MASK_VCF.smk"
include: "ng_CREATE_SNPS_TABLE.smk"
include: "ng_COMPRESS_FEAT_TABLE.smk"
include: "ng_CREATE_EXPR_TABLE.smk"

rule all:
    input:
<<<<<<< HEAD:ng_MAIN.smk
        config['tree']
'''
=======
>>>>>>> 34a99652488439fc7c9af0fa9c8ceb3cfc128396:ng_MAIN.smk
        config['expr_table'],
        config['syn_snps_table'],
        config['syn_snps_table']+'_GROUPS',
        config['syn_snps_table']+'_NON-RDNT',
        config['nonsyn_snps_table'],
        config['nonsyn_snps_table']+'_GROUPS',
        config['nonsyn_snps_table']+'_NON-RDNT'
<<<<<<< HEAD:ng_MAIN.smk
'''
=======
>>>>>>> 34a99652488439fc7c9af0fa9c8ceb3cfc128396:ng_MAIN.smk
