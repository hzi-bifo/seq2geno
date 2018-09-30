'''
rule find_best_tree:
    input:
        one_big_var_aln=TMP_D+'/OneBig.var.aln'
    output:
        tree=config['tree']
    params:
        RAXML= RAXML_EXE,
        PREFIX= 'OneBig.var',
        RESULT_D= RESULT_D,
        RAXML_OUTDIR= RESULT_D,
       	CORES= CORES
    shell:
        """
        {params[RAXML]} -T {params[CORES]} -w {params[RAXML_OUTDIR]} \
        -m GTRGAMMA -s {input} -n {params[PREFIX]} -p 1 -N 1 
        cp {params[RAXML_OUTDIR]}/RAxML_bestTree.{params[PREFIX]} {output}
        """
'''
rule find_best_tree:
    input:
        one_big_var_aln=TMP_D+'/OneBig.var.aln'
    output:
        tree=config['tree']
    params:
        FASTTREE_BIN='FastTreeMP',
        FASTTREE_OPT='-gtr -gamma -nt -quiet'
    shell:
        """
        {params.FASTTREE_BIN} {params.FASTTREE_OPT} <{input}>\
        {output.tree}
        """
    
rule postprocess_alignment:
    ## Remove the invariant sites
    input:
        one_big_aln='{TMP_D}/OneBig.aln'
    output:
        one_big_var_aln='{TMP_D}/OneBig.var.aln'
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
            "{TMP_D}/{strains}/{caller}/cons.fa", 
            TMP_D= TMP_D, strains= STRAINS, caller= 'freebayes')
    output:
        one_big_aln='{TMP_D}/OneBig.aln'

    params:
        CORES=CORES, 
        TMP_D=TMP_D+"/families", 
        STRAINS=STRAINS 
    script: 'lib/makeAlignment.py'
