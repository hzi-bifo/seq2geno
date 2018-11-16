rule for_indel_compute_gpa_raw_table:
    input:
        gffs=lambda wildcards: [os.path.join(wildcards.TMP_D, strain,
SOFTWARE['assembler'], SOFTWARE['annotator'], strain+'.gff') for strain in
DNA_READS.index.values.tolist()]
    output:
        #roary_gpa="{TMP_D}/{gene_sorter}/gene_presence_absence.csv",
        #roary_gpa_rtab='{TMP_D}/{gene_sorter}/gene_presence_absence.Rtab'
        roary_gpa="{TMP_D}/roary_for_indel/gene_presence_absence.csv",
        roary_gpa_rtab='{TMP_D}/roary_for_indel/gene_presence_absence.Rtab'
    params:
        #roary_outdir="{TMP_D}/roary",
        ROARY_BIN=software_pool.find_software('roary', target_dir=
            MAIN_VARS.library.call_lib_path()+'/roary-3.8.2/bin'),
        #ROARY_BIN='roary',
        env_cmd= ENV_POOL.activate_env_cmd('roary_env'),
        env_var_cmd= ENV_POOL.update_variables_cmd('PERL5LIB',
            MAIN_VARS.library.call_lib_path()+'/roary-3.8.2/lib/'+':$PERL5LIB'),
        cores=CORES
    shell:
        ## remove the roary folder created by snakemake first. 
        ## Otherwise, roary would create another and put all the output files in another automatically created folder
        """
        {params.env_cmd}
        {params.env_var_cmd}
        rm -r {wildcards.TMP_D}/roary_for_indel
        {params.ROARY_BIN} \
        -f {wildcards.TMP_D}/roary_for_indel -v -p \
        {params.cores} -g 100000 {input[gffs]}
        """
