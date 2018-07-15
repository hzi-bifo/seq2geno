#####
## File names are unable to determine/specify the reads files, 
## so all the file names should be carefully listed
'''
rule all:
    input:
        roary_gpa
'''
rule compute_gpa_raw_table:
    input:
        gffs=expand("{TMP_D}/prokka/{strain}.{assembler}.gff", TMP_D= TMP_D, 
            strain= STRAINS, assembler= config['assembler'])
    output:
        roary_gpa="{TMP_D}/roary/gene_presence_absence.csv"
    params:
        roary_outdir="{TMP_D}/roary",
        cores=CORES
    shell:
        ## remove the roary folder created by snakemake first. 
        ## Otherwise, roary would create another and put all the output files in another automatically created folder
        "rm -r {params.roary_outdir};"
        #"roary -f {params.roary_outdir} -e -n -v -r -p "
        "roary -f {params.roary_outdir} -v -p "
        "{params.cores} -g 100000 -z {input[gffs]};"

rule create_gff:
    input:
        assembly="{TMP_D}/{strain}/{assembler}.assem.fa"
    output:
        gff=temp("{TMP_D}/prokka/{strain}.{assembler}.gff")
    params:
        anno_outdir="{TMP_D}/prokka/{strain}",
        anno_out="{strain}.gff",
        anno_prefix='{strain}',
        cores=CORES
    shell:
        "prokka --cpus {params.cores} --force --prefix {wildcards.strain} "
        "--locustag {wildcards.strain} "
        "--outdir {params.anno_outdir} {input[assembly]};"
        "cp {params.anno_outdir}/{params.anno_out} {output[gff]}"
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
