rule import_dna_reads:
    input:
    output:
        dna_reads

rule import_dna_params:
    input:

    output:
        dna_params

rule import_reference:
    input:

    output:
        ref

rule import_rna_reads:
    input:
    
    output:
        rna_reads

rule import_rna_params:
    input:

    output:
        rna_params

rule count_gpa:
    input:
        dna_reads
        dna_params
    output:
        roary_gpa

rule detectSNPs:
    input:
        dna_reads
        dna_params
        ref
    output:
        vcf

rule compute_expressions:
    input:
        rna_reads
        rna_params
        ref
    output:
        

rule create_gpa_table:
    input:
        roary_gpa
    output:
        gpa_table

rule create_snps_table:
    input:
        vcf
    output:
        snps_table

rule infer_tree:
    input:
        vcf
        ref
    output:
        tree

rule create_expr_table:
    input:

    output:
        expr_table
             
rule make_ml_input:
    input:
        gpa_table
        snps_table
        tree
        expr_table
    output:
        ml_markdown
