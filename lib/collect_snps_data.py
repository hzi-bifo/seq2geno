import pandas as pd

def covert_snp_data(f, out_f, strains):
    # load the file 
    df=pd.read_csv(f, sep='\t', header= 0, index_col= None)

    # create feature names
    feat= ['|'.join(df.iloc[n, :6].apply(str)) for n in range(df.shape[0]) ]
    df.insert(loc=0, column= 'feature', value= feat)
    df.set_index('feature', inplace= True)

    # 
    sub_df= df[strains]
    sub_df.replace('.*\w+.*', '1', inplace= True)
    sub_df.replace('^\W$', '0', inplace= True)

    # print output
    sub_df= sub_df.transpose()
    sub_df.to_csv(out_f, sep= '\t', header= True, index= True)

if __name__=='__main__':
#    f= 'results/syn_snps.tab'
#    out_f= './test.mat'
    strains= snakemake.params['strains']
    for n in range(len(snakemake.input)):
        covert_snp_data(snakemake.input[n], snakemake.output[n], strains)
