rule cont_anc_reconstruction:
    input:
        tree_f= TREE_OUT,
        data_f= EXPR_OUT
    output: 
        output_dir=C_ANCREC_OUT
    params:
        lib_dir= os.path.join(LIB_D, 'anc_rec')
    script: 'contAncRec.R'
