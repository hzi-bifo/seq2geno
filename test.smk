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

include: "ng_DETECT_VARS.smk"
#include: "ng_CREATE_EXPR_TABLE.smk"

rule all:
    input:
       expand('{TMP_D}/{strain}/stampy/vcf.gz.tbi', TMP_D=TMP_D,
strain=STRAINS)

rule compress_feat_table:
    input: 
        F='{prefix}'
    output: 
        GROUPS='{prefix}_GROUPS',
        NONRDNT='{prefix}_NON-RDNT'
    script: 'lib/featCompress.py'
