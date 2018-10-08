rule ng_stampy_remapping:
    input:
        REF_STAMPY_INDEX=REF_FA+".stidx",
        REF_INDEXDICT=REF_FA+".sthash",
        BWA_BAM= '{TMP_D}/{strain}/stampy/pre_bwa.bam'
    output:
        STAMPY_SAM= temp('{TMP_D}/{strain}/stampy/remap.sam')
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY=STAMPY_EXE
    shell:
        """
        source activate py27
        {params[STAMPY]} \
        --readgroup=ID:{wildcards.strain},SM:{wildcards.strain}\
        -g {params.REF_PREFIX} -h {params.REF_PREFIX} \
        -t{params.CORES}  --bamkeepgoodreads -M  {input.BWA_BAM}\
        > {output.STAMPY_SAM}
        source deactivate
        """

