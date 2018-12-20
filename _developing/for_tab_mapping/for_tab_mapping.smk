rule for_tab_mapping:
    input:
        FQ=lambda wildcards: DNA_READS[wildcards.strain],
        anno_f=temp('annotations_for_snps.tab'),
        R_anno_f=temp('{TMP_D}/R_annotations.tab'),
        REF=REF_FA,
        REF_STAMPY_INDEX=REF_FA+".stidx",
        REF_STAMPY_HASH=REF_FA+".sthash",
        REF_BWA_INDEX=REF_FA+".fai",
        REF_BWA_INDEXDICT=REF_FA+".bwt"
    output:
        STAMPY_SAM='{TMP_D}/{strain}/stampy/dna_for_tab.sam'
        flc_f='{TMP_D}/{strain}/dna_for_tab.flatcount'
    params:
        prefix= lambda wildcards: os.path.join(wildcards.TMP_D, 
            wildcards.strain, 'dna_for_tab'),
        script_p='my_stampy_pipeline_PE',
        mindepth= '0'
    conda:
        ENV_FILES_POOL.find_yaml('old_mapping')
    threads: 1
    shell:
        '''     
        {params.script_p} {params.prefix} \
{params.FQ} \
{input.REF} \
{input.anno_f} \
{input.R_anno_f}
        mv {params.prefix}.sam {output.STAMPY_SAM} 
        '''
