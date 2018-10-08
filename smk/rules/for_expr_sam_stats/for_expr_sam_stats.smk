rule for_expr_sam_stats:
    input:
        STAMPY_SAM= '{TMP_D}/{strain}/stampy/rna.sam'
    output:
        STATS='{TMP_D}/{strain}/rna.stats'
    params:
        script_f= 'lib/expr/sam_statistics.pl'
    shell:
        """
        {params.script_f} -r {input.STAMPY_SAM} > {output.STATS}
        """

