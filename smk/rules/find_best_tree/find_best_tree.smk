rule find_best_tree:
    input:
        one_big_var_aln=TMP_D+'/OneBig.var.aln'
    output:
        tree=TREE_OUT
    threads: 12
    params:
        FASTTREE_BIN='FastTreeMP',
        cores= CORES,
        FASTTREE_OPT='-gtr -gamma -nt -quiet'
    shell:
        """
        export OMP_NUM_THREADS={params.cores}
        {params.FASTTREE_BIN} {params.FASTTREE_OPT} <{input}>\
        {output.tree}
        """
