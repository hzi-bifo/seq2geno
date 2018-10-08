rule postprocess_alignment:
    ## Remove the invariant sites
    input:
        one_big_aln='{TMP_D}/phylogeny/OneBig.aln'
    output:
        one_big_var_aln='{TMP_D}/phylogeny/OneBig.var.aln'
    params:
        TRIMAL_BIN='trimal',
        TRIMAL_OPT='-st 1 -gt 1 -complementary'
    shell:
        """
        {params.TRIMAL_BIN} {params.TRIMAL_OPT} \
-in {input.one_big_aln} \
-out {output.one_big_var_aln}
        """

