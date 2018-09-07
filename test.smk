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

include: "CREATE_INDEL_TABLE.smk"
include: "CREATE_SNPS_TABLE.smk"
include: "CREATE_EXPR_TABLE.salmon.smk"
include: "CREATE_GPA_TABLE.smk"
include: "COUNT_GPA.smk"
include: "CONSTRUCT_ASSEMBLY.smk"
include: "MAKE_CONS.smk"
include: "INFER_TREE.smk"
include: "DETECT_SNPS.smk"
include: "DETECT_SNPS_FOR_TABLE.smk"

rule all:
    input:
        config['expr_table']
