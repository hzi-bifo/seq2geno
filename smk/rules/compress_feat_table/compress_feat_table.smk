rule compress_feat_table:
    input: 
        F='{prefix}'
    output: 
        GROUPS='{prefix}_GROUPS',
        NONRDNT='{prefix}_NON-RDNT'
    script: 'featCompress.py'
