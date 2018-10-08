rule find_best_tree:
    input:
        one_big_var_aln=TMP_D+'/phylogeny/OneBig.var.aln'
    output:
        tree=config['tree']
    params:
        fasttree_cores= CORES,
        fasttree_bin='FastTreeMP',
        fasttree_opt='-gtr -gamma -nt -quiet'
    shell:
        """
        export OMP_NUM_THREADS={params.fasttree_cores}
        {params.fasttree_bin} {params.fasttree_opt} \
< {input.one_big_var_aln} >\
{output.tree}
        """
    
