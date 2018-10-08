rule index_bam:
    input:
        bwa_bam="{TMP_D}/{strain}/bwa/tr_bwa.bam"
    output:
        bwa_bam_index="{TMP_D}/{strain}/bwa/tr_bwa.bam.bai"
    shell:
        """
        samtools index {input.bwa_bam}
        """
