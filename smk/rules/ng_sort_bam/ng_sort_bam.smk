rule ng_sort_bam:
    input:
        paired_bam="{TMP_D}/{strain}/stampy/remap.bam"
    output:
        sorted_bam=temp("{TMP_D}/{strain}/stampy/remap_sorted.bam"),
        sorted_bam_index=temp("{TMP_D}/{strain}/stampy/remap_sorted.bam.bai")
    params:
        bamtools_bin='bamtools',
        samtools_bin='samtools'
    shell:
        """
        bamtools sort -in {input} -out {output[sorted_bam]}
        samtools index {output[sorted_bam]}
        """

