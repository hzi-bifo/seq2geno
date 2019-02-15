
def covert_snp_data(f, out_f, strains):
    import pandas as pd
    import re
    # load the file 
    df=pd.read_csv(f, sep='\t', header= 0, index_col= None)

    # create feature names
    feat= ['_'.join(df.iloc[n, :6].apply(str)) for n in range(df.shape[0]) ]
    df.insert(loc=0, column= 'feature', value= feat)
    df.set_index('feature', inplace= True)

    # 
    sub_df= df[strains]
    re_pattern= '[0-9]+\.*[0-9]*'
    sub_df= sub_df.applymap(lambda x: '1' if re.search(re_pattern, x) else '0')

    # print output
    sub_df= sub_df.transpose()
    sub_df.to_csv(out_f, sep= '\t', header= True, index= True)

if __name__=='__main__':
    strains= snakemake.params['strains']
    covert_snp_data(snakemake.input[0], snakemake.output[0], strains)
