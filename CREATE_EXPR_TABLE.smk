rule rna_mapping_bwa:
    input:
        FQ=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'rna_reads'],
        REF=REF_FA,
        REF_STAMPY_INDEX1=REF_FA+".stidx",
        REF_STAMPY_INDEX2=REF_FA+".sthash",
        REF_BWA_INDEX1=REF_FA+".bwt",
        REF_BWA_INDEX2=REF_FA+".fai"
    output:
        S_BWA_SAI= temp('{TMP_D}/{strain}/bwa/rna.sai'),
        S_BWA_BAM= temp('{TMP_D}/{strain}/bwa/rna.bam')

    params:
        BWA_OPT='-q10',
        CORES=CORES 
    shell:
        """
        bwa aln {params.BWA_OPT} -t{params.CORES} {input[REF]} {input[FQ]} \
        > {output[S_BWA_SAI]}
        bwa samse {input[REF]} {output[S_BWA_SAI]} {input[FQ]} | \
        samtools view -bS -@ {params.CORES} \
        > {output[S_BWA_BAM]}
        """

rule for_expr_stampy_remapping:
    input:
        REF_STAMPY_INDEX=REF_FA+".stidx",
        REF_INDEXDICT=REF_FA+".sthash",
        BWA_BAM= '{TMP_D}/{strain}/bwa/rna.bam'
    output:
        STAMPY_SAM= temp('{TMP_D}/{strain}/stampy/rna.sam')
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY=STAMPY_EXE
    shell:
        """
        source activate py27
        {params[STAMPY]} \
        -g {params.REF_PREFIX} -h {params.REF_PREFIX} \
        -t{params.CORES}  --bamkeepgoodreads -M  {input[BWA_BAM]}\
        > {output.STAMPY_SAM}
        source deactivate
        """

rule rna_sam2bam:
    input:
        STAMPY_SAM= '{TMP_D}/{strain}/stampy/rna.sam'
    output:
        STAMPY_BAM=temp('{TMP_D}/{strain}/stampy/rna.bam'),
        BAM=temp('{TMP_D}/{strain}/rna_sorted.bam'),
        RMDUP_SAM=temp('{TMP_D}/{strain}/rna_sorted.rmdup.sam'),
        RMDUP_BAM=temp('{TMP_D}/{strain}/rna_sorted.rmdup.bam')
    params:
        BAM_FILE_PREFIX='rna_sorted',
        CMPR_LVL=9,
        CORES=CORES
    shell:
        """
        source activate Ariane_dna
        samtools view -bS {input.STAMPY_SAM} > {output.STAMPY_BAM}
        samtools sort -@ {params.CORES} -l {params.CMPR_LVL} {output.STAMPY_BAM} \
    {wildcards.TMP_D}/{wildcards.strain}/{params.BAM_FILE_PREFIX}
        samtools rmdup -s {output.BAM} {output[RMDUP_BAM]}
        samtools view -h {output.RMDUP_BAM} > {output.RMDUP_SAM} 
        source deactivate
        """ 

'''
rule rna_mapping:
    input:
        FQ=lambda wildcards: SAMPLES_DF.loc[wildcards.sample, 'rna_reads'],
        REF=REF_FA,
        REF_STAMPY_INDEX1=REF_FA+".stidx",
        REF_STAMPY_INDEX2=REF_FA+".sthash",
        REF_BWA_INDEX1=REF_FA+".bwt",
        REF_BWA_INDEX2=REF_FA+".fai"
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY=STAMPY_EXE,
    output:
        RMDUP_SAM=temp('{TMP_D}/{strain}/rna_sorted.rmdup.sam'),
        BAM=temp('{TMP_D}/{strain}/rna_sorted.bam'),
        RMDUP_BAM=temp('{TMP_D}/{strain}/rna_sorted.rmdup.bam')
    #threads: 16
    shell:
       # 'samtools view -bS {input} |samtools sort | samtools rmdup -s |samtools view -h > {output} '
        #"source activate Ariane_dna; source activate py27;" ## needs to be dealt with wrappers
        """
        source activate py27
        {params[STAMPY]} --bwaoptions=\"-q10 {input[REF]}\" \
        -t{params[CORES]} -g {params[REF_PREFIX]} -h {params[REF_PREFIX]} \
        -M {input[FQ]} | samtools view -bS -@ {params[CORES]} | \
        samtools sort > {output[BAM]} 
        source deactivate

        source activate Ariane_dna
        samtools rmdup -s {output[BAM]} {output[RMDUP_BAM]}
        samtools view -h {output[RMDUP_BAM]} > {output[RMDUP_SAM]} 
        source deactivate
        """
'''

rule convert_to_art:
    input: 
        '{TMP_D}/{strain}/rna_sorted.rmdup.sam'
    output:
        temp('{TMP_D}/{strain}/rna_sorted.rmdup.art')
    params:
        SAM2ART='lib/sam2art.pl'
    shell:
        '{params[SAM2ART]} -s 2 {input} > {output} '

rule convert_to_sin:
    input: 
        '{TMP_D}/{strain}/rna_sorted.rmdup.sam'
    output:
        temp('{TMP_D}/{strain}/rna_sorted.rmdup.sin')
    params:
        SAM2ART='lib/sam2art.pl'
    shell:
        '{params[SAM2ART]} -s 2 -l {input} > {output} '

rule convert_to_flatcount:
    input: 
        '{TMP_D}/{strain}/rna_sorted.rmdup.sam'
    output:
        temp('{TMP_D}/{strain}/rna_sorted.rmdup.flatcount')
    params:
        SAM2ART='lib/sam2art.pl'
    shell:
        '{params[SAM2ART]} -f -s 2 {input} > {output} '

'''
rule make_annot_tab:
    input: 
        REF_GBK
    output:
        ANNO_F='Pseudomonas_aeruginosa_PA14_annotation_with_ncRNAs_07_2011_12genes.tab'
'''

rule gene_counts:
    # The required anno file is different from
    # the others in this expression pipeline
    input:
        SIN_F='{TMP_D}/{strain}/rna_sorted.rmdup.sin',
        anno_f=temp('annotations_for_snps.tab') 
    output:
        RPG_F='{TMP_D}/{strain}/rna_sorted.rmdup.rpg'
    params:
        ART2GENECOUNT='expr_original_methods/art2genecount.pl'
#        ANNO_F='Pseudomonas_aeruginosa_PA14_annotation_with_ncRNAs_07_2011_12genes.tab'
    shell:
        """
        {params.ART2GENECOUNT} -b -a {input.SIN_F} -t tab \
        -r {input.anno_f} > {output.RPG_F}
        """

rule for_expr_sam_stats:
    input:
        '{TMP_D}/{strain}/rna_sorted.rmdup.sam'
    output:
        STATS='{TMP_D}/{strain}/rna_sorted.rmdup.stats'
    params:
        script_f= 'expr_original_methods/sam_statistics.pl'
    shell:
        """
        {params.script_f} -r {input} > {output.STATS}
        """

rule for_expr_rstats:
    input:
        RPG_F='{TMP_D}/{strain}/rna_sorted.rmdup.rpg',
        STATS='{TMP_D}/{strain}/rna_sorted.rmdup.stats',
        anno_f=temp('annotations_for_expr.tab')
    output:
        RSTATS='{TMP_D}/{strain}/rna_sorted.rmdup.rstats'
    params:
        script_f='expr_original_methods/genes_statistics.R'
    shell:
        """
        {params.script_f} {input.RPG_F} {input.anno_f} {input.STATS} \
        {wildcards.TMP_D}/{wildcards.strain}/rna_sorted.rmdup
        """

rule for_expr_copy_files:
    input:
        art_f=lambda wildcards: os.path.join(TMP_D, wildcards.strain, 'rna_sorted.rmdup.art'),
        rstat_file=lambda wildcards: os.path.join(TMP_D, wildcards.strain, 'rna_sorted.rmdup.rstats')
    output:
        art_tmp_f=temp("{strain}.art"),
        rstat_tmp_f= temp("{strain}.rstats")
    shell:
        """
        ln -s {input.art_F} {output.art_tmp_f}
        ln -s {input.rstat_F} {output.rstat_tmp_f}
        """

rule for_expr_create_dict:
    input:
        art_tmp_files=expand("{strain}.art", strain=STRAINS),
        rstat_tmp_files=expand("{strain}.rstats", strain=STRAINS)
    output:
        rna_dict_f='rna_dict.txt' 
    shell:
        """
        echo {input.art_tmp_files}| \
        sed 's/\.flt\.vcf\W*/\\n/g'| \
        grep '\w' > {output.rna_dict_file}
        """

rule for_expr_art2coverage:
    input:
        art_tmp_files=expand("{strain}.art", 
            strain=STRAINS),
        rstat_tmp_files=expand("{strain}.rstats", 
            strain=STRAINS),
        rna_dict_f='rna_dict.txt' 
    output:
        '{TMP_D}/coverage_cut1.txt'
    params:
        script_f='expr_original_methods/art2cov.py'
    shell:
        """
        {params.script_f} \
        -f {input.rna_dict_f} -s SE -c 1 \
        -o {output}
        """

rule write_dict_file:
    input:
        RPG_FILES=expand('{TMP_D}/{strain}/rna_sorted.rmdup.rpg',
strain=STRAINS, TMP_D= TMP_D)
    output:
        rpg_dict_f= TMP_D+'/rpg_dict'
    params:
        tmp_d= TMP_D,
        strains= STRAINS,
        target_filename='rna_sorted.rmdup.rpg'
    shell:
        """
        echo 'name  path' > {output.rpg_dict_f}
        ls {params.tmp_d}/*/{params.target_filename}|\
        awk -v pwd=$PWD \
        -F"/" '{{print $2"\t"pwd"/"$0}}' >> {output.rpg_dict_f}
        """

rule for_expr_create_annot:
    input:
        ref_gbk=config['reference_annotation']
    output:
        anno_f=temp('annotations_for_expr.tab')
    script:'lib/create_dict_for_expr.py'

rule create_and_make_expr_table:
    input:
        rpg_dict_f= TMP_D+'/rpg_dict',
        anno_f=temp('annotations_for_expr.tab')
    output:
        expr_table=config['expr_table']
    script: 'lib/collect_rpg_data.R'
