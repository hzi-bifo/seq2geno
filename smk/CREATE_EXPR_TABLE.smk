rule create_and_make_expr_table:
    input:
        RPG_FILES=expand('{TMP_D}/{strain}/rna.rpg', strain=STRAINS, TMP_D=
TMP_D)
    output:
        expr_table=config['expr_table']
    params:
        tmp_d= TMP_D,
        strains= STRAINS
    run:
        import pandas as pd
        import re

        def obtain_genecount_series(strain, rpg_f_with_feat):
            rpg_df=pd.read_csv(rpg_f_with_feat, sep='\t', header= None,
index_col= None)
            colnames= ['feat', 'gene_count', 'antisensecount']
            rpg_df= rpg_df.rename(columns= {n: colnames[n] for n in
range(len(colnames))}).set_index('feat')
            return(rpg_df['gene_count'].rename(strain))
        # the rpg files
        rpg_files= {strain: os.path.join(params['tmp_d'], strain, 'rna.rpg') for strain in params['strains']}
        # integrate into one table
        rpg_df= pd.concat([obtain_genecount_series(s, rpg_files[s]) for s in
rpg_files],axis= 1).T
        rpg_df.to_csv(output['expr_table'], sep= '\t', header= True, index= True)

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

rule for_expr_create_annot:
    input:
        ref_gbk=config['reference_annotation']
    output:
        R_anno_f=temp('{TMP_D}/R_annotations.tab')
    script:'lib/create_R_anno.py'

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


rule convert_to_art:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/stampy/rna.sam'
    output:
        ART=temp('{TMP_D}/{strain}/rna.art')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params.SAM2ART} -s 2 {input.STAMPY_SAM} > {output} '

rule convert_to_sin:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/stampy/rna.sam'
    output:
        SIN=temp('{TMP_D}/{strain}/rna.sin')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params.SAM2ART} -s 2 -l {input.STAMPY_SAM}  > {output} '

rule convert_to_flatcount:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/stampy/rna.sam'
    output:
        FLTCNT=temp('{TMP_D}/{strain}/rna.flatcount')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params.SAM2ART} -f -s 2 {input.STAMPY_SAM}  > {output} '


rule rna_mapping:
    input:
        FQ=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'rna_reads'],
        REF=REF_FA,
        REF_STAMPY_INDEX1=REF_FA+".stidx",
        REF_STAMPY_INDEX2=REF_FA+".sthash",
        REF_BWA_INDEX1=REF_FA+".bwt",
        REF_BWA_INDEX2=REF_FA+".fai"
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY=STAMPY_EXE
    output:
        STAMPY_SAM= temp('{TMP_D}/{strain}/stampy/rna.sam')
    #threads: 16
    shell:
       # 'samtools view -bS {input} |samtools sort | samtools rmdup -s |samtools view -h > {output} '
        #"source activate Ariane_dna; source activate py27;" ## needs to be dealt with wrappers
        """
        source activate Ariane_dna
        {params[STAMPY]} --bwaoptions=\"-q10 {input[REF]}\" \
-g {params.REF_PREFIX} \
-h {params.REF_PREFIX} \
-M {input.FQ} > {output.STAMPY_SAM} 
        source deactivate
        """

