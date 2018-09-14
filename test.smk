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

include: "ng_MASK_VARS.smk"
include: "ng_MERGE_VARS.smk"
include: "ng_DETECT_VARS.smk"
#include: "ng_CREATE_EXPR_TABLE.smk"

rule all:
    input:
        os.path.join(TMP_D, 'stampy', 'multisample_vcf.gz')
#        coding_vcf_gz=os.path.join(TMP_D, "multi_sample_vcf.coding.gz"),
#        igr_vcf_gz=os.path.join(TMP_D, "multi_sample_vcf.igr.gz")

rule compress_feat_table:
    input: 
        F='{prefix}'
    output: 
        GROUPS='{prefix}_GROUPS',
        NONRDNT='{prefix}_NON-RDNT'
    script: 'lib/featCompress.py'
