'''
Compress binary features by their patterns among the samples.
A representiative is selected for each pattern. The new feature name is composed of the representative feature 
and the size of group.
'''

import pandas as pd
import sys

def check_input(df):
    ## check if the features (only 1/0) are binary and contain no NA values
    check_df= df.fillna('NA')
    outcome= True if ((check_df != '1') & (check_df != '0')).sum().sum() == 0 else False
    if not outcome:
        sys.exit("File unrecognizable: only '1' and '0' are acceptable and shouldn't contain NA values")

f='results/annot.tab'
group_out_f= 'results/annot.tab'+'_GROUPS'
comp_out_f= 'results/annot.tab'+'_NONRDNT'

## read table
df= pd.read_csv(f, sep= '\t', header= 0, index_col= 0, dtype= 'str')
check_input(df)

## group the features
t_df= df.transpose() # features in rows
grouped_t_df= t_df.groupby(t_df.columns.values.tolist())

## list the groups
## and convert the feature names
new_feat_dict= {}
group_dict= {}
for name, group in grouped_t_df:
    grouped_t_df.get_group(name)
    repre= str(group.index.tolist()[0])
    new_repre= '|'.join([repre, str(len(group))]) 
    feat= t_df.loc[repre,:]
    
    new_feat_dict[new_repre]= feat
    group_dict[new_repre]= group.index.tolist()

## print the compressed features
new_df= pd.DataFrame(data= new_feat_dict) # features in columns
new_df.to_csv(comp_out_f, sep= '\t', header= True, 
    index= True, index_label= '', encoding= 'utf-8')

## print the cluster information
with open(group_out_f, 'w') as group_out_fh:
    for new_repre in group_dict:
        group_out_fh.write('{}:\t{}\n'.
	format(new_repre, '\t'.join(group_dict[new_repre])))

