rule gene_counts:
    # The required anno file is different from
    # the others in this expression pipeline
    input:
        SIN='{TMP_D}/{strain}/rna.sin',
        anno_f='annotations_for_snps.tab' 
    output:
        RPG_F='{TMP_D}/{strain}/rna.rpg'
    params:
        ART2GENECOUNT='lib/expr/art2genecount.pl'
    shell:
        """
        {params.ART2GENECOUNT} \
-a {input.SIN} -t tab \
-r {input.anno_f} > {output.RPG_F}
        """


