rule find_best_tree:
    input:
        one_big_var_aln=TMP_D+'/OneBig.var.aln'
    output:
        tree=TREE_OUT
    params:
        FASTTREE_BIN='FastTreeMP',
        FASTTREE_OPT='-gtr -gamma -nt -quiet'
    shell:
        """
        {params.FASTTREE_BIN} {params.FASTTREE_OPT} <{input}>\
        {output.tree}
        """
