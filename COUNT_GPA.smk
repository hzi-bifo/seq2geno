#####
## File names are unable to determine/specify the reads files, 
## so all the file names should be carefully listed

rule all:
    input:
        roary_gpa

rule compute_gpa_raw_table:
    input:
        gffs=[...]
    output:
        roary_gpa
    params:
        roary_dir
        cores=
    shell:
        ## remove the roary folder created by snakemake first. Otherwise, roary would create another and put all the output files there
        "rm -r {params.roary_dir};"
        "roary -f {params.roary_dir} -e -n -v -r -p {params.cores} -g 100000 -z {input[gffs]};"

rule create_gff:
    input:
        assembly=
    output:
        gff=
    params:
        anno_dir=
        cores=
    shell:
        "prokka --cpus {params.cores} --force --prefix {wildcards.strain} --outdir {params.anno_dir} {input[assembly]}"

rule create_assembly:
    input:
        reads1=
        reads2=
    output:
        assembly=
    params:
        assem_dir
        cores
    shell:
        "spades.py --careful --threads {params.cores} -o {params.assem_dir}/assemblies/{wildcards.strain} -1 {input[reads1]} -2 {input[reads2]}"

rule load_reads:
    output:
        reads1=
        reads2=
'''
#####
rule convert_gpa_table:
    input:
        TMP_D+'/roary/gene_presence_absence.csv'
    output:
        'gpa.mat'
    script: 'roary_gpa2bin.R'

rule run_roary:
    input:
        expand(TMP_D+"/annotations/{strain}/{strain}.gff", strain= STRAINS)
    output:
        TMP_D+'/roary/gene_presence_absence.csv'
    shell:
        "rm -r {TMP_D}/roary"
        "roary -f {TMP_D}/roary -e -n -v -r -p CORES -g 100000 -z {input}"

rule run_prokka:
    input:
        scaf_f=TMP_D+'/assemblies/{strain}/scaffolds.fasta'
    output:
        TMP_D+"/annotations/{strain}/{strain}.gff"
    params:
        OUT_DIR= TMP_D+"/annotations/{strain}"
    shell:
        "prokka --cpus {CORES} --force --prefix {wildcards.strain} --outdir {params.OUT_DIR} {input}"

rule make_assembly:
    input:
#        lambda wildcards: [SAMPLES_DF.loc[wildcards.strain, 'reads1'], SAMPLES_DF.loc[wildcards.strain, 'reads2']]
        reads1
        reads2
    output:
        scaf_f=TMP_D+'/assemblies/{strain}/scaffolds.fasta'
    shell:
        "spades.py --careful --threads {CORES} -o {TMP_D}/assemblies/{wildcards.strain} -1 {input[reads1]} -2 {input[reads2]}"
'''
