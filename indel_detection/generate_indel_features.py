import pandas as pd
# list all the indel detection outcomes
list_f=snakemake.input['indel_list']
# create the indel table
files= [l.strip() for l in open(list_f, 'r')]
indel_df=pd.concat([pd.read_csv(f, sep= '\t', index_col= 0) for f in files],
        axis=1)
indel_df.rename(index=str, inplace=True, columns={g: g+'_indel' for g in
    indel_df.columns.values.tolist()})
# count the frequency
# the 0s are presence
freq=indel_df.shape[1]-indel_df.sum(axis=1)
# write the files
indel_df.to_csv(snakemake.output['annot_f'], sep= '\t')
freq.to_csv(snakemake.output['indel_stat_f'], sep= '\t')
