
rule for_tab_stampy_mapping:
    input:
        FQ=lambda wildcards: DNA_READS[wildcards.strain],
        REF=REF_FA,
        REF_STAMPY_INDEX=REF_FA+".stidx",
        REF_STAMPY_HASH=REF_FA+".sthash",
        REF_BWA_INDEX=REF_FA+".fai",
        REF_BWA_INDEXDICT=REF_FA+".bwt"
    output:
        STAMPY_SAM='{TMP_D}/{strain}/stampy/dna_for_tab.sam'
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY='stampy.py'
        #STAMPY=STAMPY_EXE
    conda:
        ENV_FILES_POOL.find_yaml('old_mapping')
    shell:
        '''
        {params.STAMPY} --bwaoptions=\"-q10 {input.REF}\" \
-g {params.REF_PREFIX} \
-h {params.REF_PREFIX} \
-M {input.FQ} > {output.STAMPY_SAM} 
        '''    
