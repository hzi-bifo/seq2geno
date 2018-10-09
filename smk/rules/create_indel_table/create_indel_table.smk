rule create_indel_table:
    input:
        indel_list= os.path.join(TMP_D, 'indel.list')
    output:
        annot_f= INDEL_OUT,
        indel_stat_f='indels_stats.txt'
    script: 'generate_indel_features.py'
