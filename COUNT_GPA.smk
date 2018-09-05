#####
## File names are unable to determine/specify the reads files, 
## so all the file names should be carefully listed
rule compute_gpa_raw_table:
    input:
        gffs=lambda wildcards: [os.path.join(wildcards.TMP_D, strain, SOFTWARE['assembler'], SOFTWARE['annotator'], strain+'.gff') for strain in STRAINS]
    output:
        roary_gpa="{TMP_D}/{gene_sorter}/gene_presence_absence.csv",
        roary_gpa_rtab='{TMP_D}/{gene_sorter}/gene_presence_absence.Rtab'
    params:
        #roary_outdir="{TMP_D}/roary",
        ROARY_BIN='roary',
        cores=CORES
    shell:
        ## remove the roary folder created by snakemake first. 
        ## Otherwise, roary would create another and put all the output files in another automatically created folder
        """

        rm -r {wildcards.TMP_D}/{wildcards.gene_sorter}
        {params.ROARY_BIN} -b \"blastp -qcov_hsp_perc 95 \" \
        -f {wildcards.TMP_D}/{wildcards.gene_sorter} -v -p \
        {params.cores} -g 100000 {input[gffs]}

        """

rule duplicate_gffs:
    ## merely for roary
    input:
        gff_output="{TMP_D}/{strain}/{assembler}/{annotator}/de_novo.gff"
    output:
        gff_output_copy=temp("{TMP_D}/{strain}/{assembler}/{annotator}/{strain}.gff")
    shell:
        """
        cp {input[gff_output]} {output[gff_output_copy]}
        """
        
rule create_gff:
    input:
        assembly="{TMP_D}/{strain}/{assembler}/scaffolds.fasta"
    output:
        ffn_output="{TMP_D}/{strain}/{assembler}/{annotator}/de_novo.ffn",## required by the indel table
        gff_output="{TMP_D}/{strain}/{assembler}/{annotator}/de_novo.gff"
    params:
        #anno_outdir="{TMP_D}/prokka/{strain}",
        #anno_out="{strain}.gff",
        anno_prefix='de_novo',
        #anno_prefix='{strain}',
        cores=CORES
    shell:
        """

        prokka --cpus {params.cores} --force --prefix {params.anno_prefix} \
        --locustag {wildcards.strain} \
        --outdir {wildcards.TMP_D}/{wildcards.strain}/{wildcards.assembler}/{wildcards.annotator} {input[assembly]}

        """
#        "cp {params.anno_outdir}/{params.anno_out} {output[gff]}"
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
