rule dif_xpr_analysis:
    input:
        phe_table= PHE_TABLE_F,
        expr_table= EXPR_OUT
    output:
        dif_xpr_out=DIF_XPR_OUT 
    params:
        lib_dir= os.path.join(LIB_D, 'dif_xpr'),
        alpha_cutoff=DIFXPR_ALPHA,
        lfc_cutoff= DIFXPR_LFC, 
        cores= CORES, 
    script: 'diffEpr_analysis.R'
        
