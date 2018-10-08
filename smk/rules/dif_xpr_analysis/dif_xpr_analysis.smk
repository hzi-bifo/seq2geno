rule dif_xpr_analysis:
    input:
        phe_table=config['phe_table'],
        expr_table= config['expr_table']
    output:
        dif_xpr_out= directory(os.path.join(TMP_D, 'deseq2'))
    params:
        lib_dir= os.path.join(LIB_D, 'dif_xpr'),
        alpha_cutoff=config['dif_xpr_alpha'],#from config
        lfc_cutoff= config['dif_xpr_lfc'], #from config
        cores= CORES, 
    script: 'diffEpr_analysis.R'
        
