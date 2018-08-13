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
RESULT_D=config['result_d']
#SOFTWARE={'mapper': 'bwa'}
SOFTWARE= config['software']
SOFTWARE['annotator']= 'prokka'
SOFTWARE['gene_sorter']= 'roary'
print(SAMPLES_DF)
print(SOFTWARE)

rule all:
    input:
        FAM_INDELS_TXT=os.path.join('.','test_indels.txt'),

rule vcf_to_indels_per_fam:
    input:
        FAM_VCF=os.path.join('.', '{fam}.vcf')       
    output:
        FAM_INDELS_TXT=os.path.join('.', '{fam}_indels.txt'),
    params:
        cores= CORES,
       # fam=wildcards.fam,
       # vcf2indel_script='indel_detection/vcf2indel.rewrite.py',
        working_dir='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3',
        strains_perc_cutoff= 0.5,
        len_cutoff= 8
    script: 'indel_detection/vcf2indel.rewrite.py'
        
