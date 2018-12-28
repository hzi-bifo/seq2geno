rule compute_gpa_table:
    input:
        roary_gpa=TMP_D+"/roary/gene_presence_absence.csv"
    output:
        gpa_table=GPA_OUT
    params: 
        strains= DNA_READS.index.values.tolist()
    script: 'roary_gpa2bin.py'
