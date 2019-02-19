import pandas as pd
import textwrap

strains=snakemake.params['strains']
rename_dict_f=snakemake.params['rename_dict']
gpa_f= snakemake.input['roary_gpa']
out_f=snakemake.output['gpa_table']

## read roary result and identify single-copy genes
df=pd.read_csv(gpa_f, sep= ',',
    header= 0, index_col= 0, quotechar= '"', low_memory=False)

## read the dictionary
rename_d= {}
with open(rename_dict_f, 'r') as rename_dict_fh:
    for l in rename_dict_fh:
        d= l.strip('\n').split('\t')
        rename_d[d[0]]= d[1]

## filter and convert the states
sub_df=df.loc[:, strains]
sub_df= sub_df.applymap(lambda x: '1' if (not pd.isna(x)) else '0')
sub_df= sub_df.transpose() # strains in rows
sub_df.rename(columns=rename_d, inplace= True)

## print 
sub_df.to_csv(out_f, sep= '\t', index= True, index_label='')
