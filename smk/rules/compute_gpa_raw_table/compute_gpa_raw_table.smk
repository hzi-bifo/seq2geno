rule compute_gpa_raw_table:
    input:
        gffs=lambda wildcards: [os.path.join(wildcards.TMP_D, strain,
SOFTWARE['assembler'], SOFTWARE['annotator'], strain+'.gff') for strain in
DNA_READS.index.values.tolist()]
    output:
        roary_gpa="{TMP_D}/roary/gene_presence_absence.csv",
        roary_gpa_rtab='{TMP_D}/roary/gene_presence_absence.Rtab'
    params:
        #roary_outdir="{TMP_D}/roary",
        #ROARY_BIN=software_pool.find_software('roary', include_env= True),
        ROARY_BIN='roary',
    threads:1
    conda: ENV_FILES_POOL.find_yaml('roary_env')
    shell:
        ## remove the roary folder created by snakemake first. 
        ## Otherwise, roary would create another and put all the output files in another automatically created folder
        '''
        echo $PERL5LIB
        rm -r {wildcards.TMP_D}/roary
        {params.ROARY_BIN} \
-f {wildcards.TMP_D}/roary -v -p \
{threads} -g 100000 {input.gffs}
        '''
