rule ng_create_consensus_sequence:
    ## Create consensus sequences of the coding regions 
    input:
        ref_target_seqs="{TMP_D}/reference.target_regions.fa",
        coding_vcf_gz="{TMP_D}/freebayes/multisample.vcf.coding.gz",
        coding_vcf_gz_index="{TMP_D}/freebayes/multisample.vcf.coding.gz.tbi"
    output:
        cons_coding_seqs="{TMP_D}/{strain}/cons.fa"
    shell:
        """
        bcftools consensus --sample {wildcards.strain} \
-f {input.ref_target_seqs} \
        {input.coding_vcf_gz} > {output.cons_coding_seqs}
        """

