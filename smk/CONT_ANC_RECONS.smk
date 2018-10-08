rule cont_anc_reconstruction:
    input:
        tree_f= config['tree'],
        data_f= config['expr_table']
    output: 
        output_dir=directory(os.path.join(TMP_D, 'fastAnc'))
    params:
        lib_dir= os.path.join(LIB_D, 'anc_rec')
    script: '../lib/anc_rec/contAncRec.R'
