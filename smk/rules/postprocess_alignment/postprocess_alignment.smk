rule postprocess_alignment:
    ## Remove the invariant sites
    input:
        one_big_aln='{TMP_D}/OneBig.aln'
    output:
        one_big_var_aln='{TMP_D}/OneBig.var.aln'
    params:
        TRIMAL_BIN='trimal',
        TRIMAL_OPT='-st 1 -gt 1 -complementary'
    threads: 1
    shell:
        """
        {params.TRIMAL_BIN} {params.TRIMAL_OPT} \
-in {input.one_big_aln} \
-out {output.one_big_var_aln}
        """

