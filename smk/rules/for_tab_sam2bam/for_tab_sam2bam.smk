rule for_tab_sam2bam:
    input:
        STAMPY_SAM= temp('{TMP_D}/{strain}/stampy/dna_for_tab.sam')
    output:
        bam= temp('{TMP_D}/{strain}/stampy/tab_dna.bam'),
        bam_index=temp("{TMP_D}/{strain}/stampy/tab_dna.bam.bai")
    params: 
        output_prefix= lambda wildcards: os.path.join(wildcards.TMP_D,
wildcards.strain, 'stampy', 'tab_dna')
    shell:
        """
        source activate Ariane_dna
        samtools view -bS \
{input.STAMPY_SAM} > {output.bam}
        samtools sort {output.bam} \
{params.output_prefix}
        samtools index {output.bam}
        source deactivate
        """

