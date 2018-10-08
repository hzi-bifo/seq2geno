rule create_consensus_sequence:
    ## Create consensus sequences of the coding regions 
    input:
        ref_target_seqs="{TMP_D}/reference.target_regions.fa",
        vcf_gz="{TMP_D}/{strain}/{mapper}/vcf.gz",
        vcf_gz_index="{TMP_D}/{strain}/{mapper}/vcf.gz.tbi"
    output:
        cons_coding_seqs="{TMP_D}/{strain}/{mapper}/cons.fa"
    shell:
        """
        bcftools consensus -f {input[ref_target_seqs]} \
        {input[vcf_gz]} > {output}
        """

