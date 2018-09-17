'''
Purpose:
Remove the redundant features by binary patterns among the samples

Output:
        GROUPS: The clustering result, where each row includes those sharing
the same pattern
        NONRDNT: The table with non-redundant features set
'''
rule compress_feat_table:
    input: 
        F='{prefix}'
    output: 
        GROUPS='{prefix}_GROUPS',
        NONRDNT='{prefix}_NON-RDNT'
    script: 'lib/featCompress.py'
