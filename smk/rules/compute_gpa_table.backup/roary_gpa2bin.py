import pandas as pd
import textwrap

strains=snakemake.params['strains']
gpa_f= snakemake.input['roary_gpa']
out_f=snakemake.output['gpa_table']

## read roary result and identify single-copy genes
df=pd.read_csv(gpa_f, sep= ',',
        header= 0, index_col= 0, quotechar= '"', low_memory=False)

## filter and convert the states
sub_df=df.loc[:, strains]
sub_df= sub_df.applymap(lambda x: '1' if (not pd.isna(x)) else '0')
sub_df= sub_df.transpose() # strains in rows

## print 
sub_df.to_csv(out_f, sep= '\t', index= True, index_label='')
