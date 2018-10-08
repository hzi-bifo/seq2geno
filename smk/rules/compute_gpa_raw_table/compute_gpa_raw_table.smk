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

