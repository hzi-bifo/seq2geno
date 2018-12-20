rule for_snps_convert_to_flatcount:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/stampy/dna_for_tab.sam'
    output:
        FLTCNT='{TMP_D}/{strain}/dna_for_tab.flatcount'
    params:
        script_f='sam2art.pl'
    conda:
        ENV_FILES_POOL.find_yaml('old_mapping')
    shell:
        '{params.script_f} -f -s 2 -p {input.STAMPY_SAM}  > {output} '

