import pandas as pd

#def covert_snp_data(f, out_f, strains):
def covert_snp_data(f, out_f):
    # load the file 
    df=pd.read_csv(f, sep='\t', header= 0, index_col= None)

    # create feature names
    feat= ['|'.join(df.iloc[n, :6].apply(str)) for n in range(df.shape[0]) ]
    df.insert(loc=0, column= 'feature', value= feat)
    df.set_index('feature', inplace= True)

    # 
    sub_df=df.iloc[:, 6:]
    strains=sub_df.columns.values.tolist()
    sub_df.replace('.*\w+.*', '1', inplace= True)
    sub_df.replace('^\W$', '0', inplace= True)

    # print output
    sub_df= sub_df.transpose()
    sub_df.to_csv(out_f, sep= '\t', header= True, index= True)

if __name__=='__main__':
    f= 'tmp/non-syn_SNPs.tab'
    out_f= './test.mat'
    covert_snp_data(f, out_f)
