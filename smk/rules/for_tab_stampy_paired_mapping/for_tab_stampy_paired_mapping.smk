rule for_tab_stampy_paired_mapping:
    input:
        FQ1=lambda wildcards: DNA_READS[wildcards.strain][0],
        FQ2=lambda wildcards: DNA_READS[wildcards.strain][1],
        REF=REF_FA,
        REF_BWA_INDEX=REF_FA+".fai",
        REF_BWA_INDEXDICT=REF_FA+".bwt"
    output:
        STAMPY_RAW_SAM= temp('{TMP_D}/{strain}/stampy/dna_for_tab.2.sam')
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY=STAMPY_EXE
    shell:
        """
        source activate old_mapping
        {params.STAMPY} --bwaoptions=\"-q10 {input.REF}\" \
-g {params.REF_PREFIX} \
-h {params.REF_PREFIX} \
-M {input.FQ1} {input.FQ2} > {output.STAMPY_RAW_SAM} 
        source deactivate
        """
