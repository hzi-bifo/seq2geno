rule for_expr_rstats:
    input:
        RPG_F='{TMP_D}/{strain}/rna.rpg',
        STATS='{TMP_D}/{strain}/rna.stats',
        R_anno_f='{TMP_D}/R_annotations.tab'
    output:
        RSTATS='{TMP_D}/{strain}/rna'
    params:
        script_f='lib/expr/genes_statistics.R'
    shell:
        """
        {params.script_f} {input.RPG_F} {input.R_anno_f} {input.STATS} \
        {wildcards.TMP_D}/{wildcards.strain}/rna
        """

