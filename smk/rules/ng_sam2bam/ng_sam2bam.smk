rule ng_sam2bam:
    input:
        STAMPY_SAM= temp('{TMP_D}/{strain}/stampy/remap.sam')
    output:
        STAMPY_BAM= temp('{TMP_D}/{strain}/stampy/remap.bam')
    params:
        CORES=CORES
    shell:
        """
        source activate old_mapping
        samtools view -bS -@ {params.CORES} \
        {input.STAMPY_SAM} > {output.STAMPY_BAM}
        source deactivate
        """

