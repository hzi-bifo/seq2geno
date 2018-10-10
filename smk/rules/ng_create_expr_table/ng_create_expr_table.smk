rule ng_create_expr_table:
    input :
        SALMON_OUTPUTS= expand(directory('{TMP_D}/{strain}/salmon'),
TMP_D=TMP_D, strain= RNA_READS.index.values.tolist())
    output:
        expr_table=EXPR_OUT
    params:
        strains= RNA_READS.index.values.tolist(),
        tmp_d= TMP_D,
        salmon_raw_output= 'quant.sf'
    run:
        import pandas as pd
        import os

        files_dict= {s: os.path.join(params.tmp_d, s, 'salmon', 'quant.sf') for
s in params.strains}
        ## open and read the Salmon outputs
        rnum_dict= {s: pd.read_csv(files_dict[s], sep= '\t', header= 0,
            index_col= 0, dtype= str)['NumReads'] for s in params.strains}
        ## convert to a matrix
        rnum_df= pd.DataFrame(data= rnum_dict).transpose() # strains in rows
        rnum_df.to_csv(output.expr_table,
            sep= '\t',
            na_rep= 'NA',
            header=True, index= True)

