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

rule create_coding_regions_aln:
    ## Sort sequences by gene family, align each family, and concatenate the
    ## consensus sequence alignments
    input:
        cons_coding_seqs_every_strain=expand(
            "{TMP_D}/{strain}/cons.fa", 
            TMP_D= TMP_D, strain= STRAINS)
    output:
        one_big_aln='{TMP_D}/phylogeny/OneBig.aln'
    params:
        CORES=CORES, 
        TMP_D=TMP_D+"/phylogeny/families", 
        STRAINS=STRAINS 
    script: 'lib/makeAlignment.py'
