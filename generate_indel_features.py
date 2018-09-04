import pandas as pd

# list all the indel detection outcomes
list_f='seq2geno_temp/indel.list'
out_annot_f= 'results/annot.tab'
out_stat_f= 'indels_stats.txt'
# create the indel table
files= [l.strip() for l in open(list_f, 'r')]
indel_df=pd.concat([pd.read_csv(f, sep= '\t', index_col= 0, dtype=str) for f in files],
        axis=1)
# those without indels were not listed in the indel list
indel_df.fillna('0', inplace= True)
indel_df.rename(index=str, inplace=True, columns={g: g+'_indel' for g in
    indel_df.columns.values.tolist()})
# count the frequency
# the 0s are presence
freq=indel_df.shape[1]-indel_df.sum(axis=1)
# write the files
indel_df.to_csv(out_annot_f, sep= '\t')
freq.to_csv(out_stat_f, sep= '\t')
