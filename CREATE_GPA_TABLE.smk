rule all:
    input:
        gpa_table

rule compute_gpa_table:
    input:
        ## how to define/find the file
        roary_gpa
    output:
        gpa_table
    script: 'roary_gpa2bin.R'

rule load_roary
    output:
        roary_gpa

