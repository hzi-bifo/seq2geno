rule compute_gpa_table:
    input:
        ## how to define/find the file
        roary_gpa=TMP_D+"/roary/gene_presence_absence.csv"
    output:
        gpa_table=config['gpa_table']
    script: 'roary_gpa2bin.R'
