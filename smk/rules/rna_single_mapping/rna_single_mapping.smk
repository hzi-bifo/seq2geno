rule rna_single_mapping:
    input:
        FQ=lambda wildcards: RNA_READS[wildcards.strain][0],
        REF=REF_FA,
        REF_STAMPY_INDEX1=REF_FA+".stidx",
        REF_STAMPY_INDEX2=REF_FA+".sthash",
        REF_BWA_INDEX1=REF_FA+".bwt",
        REF_BWA_INDEX2=REF_FA+".fai"
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY=STAMPY_EXE
    output:
        STAMPY_RAW_SAM= temp('{TMP_D}/{strain}/stampy/rna.1.sam')
    shell:
        """
        source activate old_mapping
        {params[STAMPY]} --bwaoptions=\"-q10 {input[REF]}\" \
-g {params.REF_PREFIX} \
-h {params.REF_PREFIX} \
-M {input.FQ} > {output.STAMPY_RAW_SAM} 
        source deactivate
        """

