rule create_expr_table:
    input:
        RPG_FILES=expand('{TMP_D}/{strain}/rna.rpg', 
        strain=RNA_READS.index.values.tolist(), 
        TMP_D=TMP_D)
    output:
        expr_table=EXPR_OUT
    params:
        tmp_d= TMP_D,
        strains= RNA_READS.index.values.tolist()
    run:
        import pandas as pd
        import re

        def obtain_genecount_series(strain, rpg_f_with_feat):
            rpg_df=pd.read_csv(rpg_f_with_feat, sep='\t', header= None,
index_col= None)
            colnames= ['feat', 'gene_count', 'antisensecount']
            rpg_df= rpg_df.rename(columns= {n: colnames[n] for n in
range(len(colnames))}).set_index('feat')
            return(rpg_df['gene_count'].rename(strain))
        # the rpg files
        rpg_files= {strain: os.path.join(params['tmp_d'], strain, 'rna.rpg') for strain in params['strains']}
        # integrate into one table
        rpg_df= pd.concat([obtain_genecount_series(s, rpg_files[s]) for s in
rpg_files],axis= 1).T
        rpg_df.to_csv(output['expr_table'], sep= '\t', header= True, index= True)

