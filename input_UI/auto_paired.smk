import pandas as pd
import os
import re
def parse_samples(f):
    samples_df= pd.read_csv(f, sep= '\t', header=0, dtype=
str).set_index("strain",drop=False)
    print(samples_df)
    if samples_df.isna().any(axis=None) :
        exit('ERROR: NA values found in the samples table')
        
    samples_df['rna_reads']= samples_df['rna_reads'].apply(lambda x:
        re.split('\s*,\s*', x))
    samples_df['dna_reads']= samples_df['dna_reads'].apply(lambda x:
        re.split('\s*,\s*', x))
    
    ## detect the number of reads files
    samples_df['rna_reads_layout']= samples_df['rna_reads'].apply(lambda x: len(x))
    samples_df['dna_reads_layout']= samples_df['dna_reads'].apply(lambda x: len(x))

    ## validate
    if ((samples_df['rna_reads_layout'].apply(lambda x: not (x in [1,2])).any())
| (samples_df['dna_reads_layout'].apply(lambda x: not (x in [1,2])).any())):
        exit('ERROR: Unusual number of reads detected')

    return(samples_df)
    
configfile: "config.yaml"
SAMPLES_DF=parse_samples(config["samples"])
STRAINS=SAMPLES_DF['strain'].tolist()

print(SAMPLES_DF)
