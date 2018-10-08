rule rna_paired_mapping:
    input:
        FQ1=lambda wildcards: RNA_READS[wildcards.strain][0],
        FQ2=lambda wildcards: RNA_READS[wildcards.strain][1],
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
        STAMPY_RAW_SAM= temp('{TMP_D}/{strain}/stampy/rna.2.sam')
    shell:
        """
        source activate Ariane_dna
        {params[STAMPY]} --bwaoptions=\"-q10 {input[REF]}\" \
-g {params.REF_PREFIX} \
-h {params.REF_PREFIX} \
-M {input.FQ1} {input.FQ2} > {output.STAMPY_RAW_SAM} 
        source deactivate
        """

