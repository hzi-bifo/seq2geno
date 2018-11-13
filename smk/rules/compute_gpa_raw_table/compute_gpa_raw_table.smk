rule compute_gpa_raw_table:
    input:
        gffs=lambda wildcards: [os.path.join(wildcards.TMP_D, strain,
SOFTWARE['assembler'], SOFTWARE['annotator'], strain+'.gff') for strain in
DNA_READS.index.values.tolist()]
    output:
        roary_gpa="{TMP_D}/{gene_sorter}/gene_presence_absence.csv",
        roary_gpa_rtab='{TMP_D}/{gene_sorter}/gene_presence_absence.Rtab'
    params:
        #roary_outdir="{TMP_D}/roary",
        #ROARY_BIN='roary',
        ROARY_BIN=software_pool.find_software('roary', include_env= True),
        cores=CORES
    shell:
        ## remove the roary folder created by snakemake first. 
        ## Otherwise, roary would create another and put all the output files in another automatically created folder
        """
        source activate roary_env
        export PERL5LIB=$CONDA_PREFIX/lib/perl5/site_perl/5.22.0/:$PERL5LIB
        rm -r {wildcards.TMP_D}/{wildcards.gene_sorter}
        {params.ROARY_BIN} \
        -f {wildcards.TMP_D}/{wildcards.gene_sorter} -v -p \
        {params.cores} -g 100000 {input[gffs]}
        """

