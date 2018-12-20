rule for_tab_sam2fltvcf:
    input:
        REF=REF_FA,
        REF_STAMPY_INDEX=REF_FA+".stidx",
        REF_STAMPY_HASH=REF_FA+".sthash",
        REF_BWA_INDEX=REF_FA+".fai",
        REF_BWA_INDEXDICT=REF_FA+".bwt",
        STAMPY_SAM='{TMP_D}/{strain}/stampy/dna_for_tab.sam'
    output:
        snp_vcf_file='{TMP_D}/{strain}/samtools/tab_dna.flt.vcf'
    params:
        prefix= lambda wildcards: os.path.join(wildcards.TMP_D, 
            wildcards.strain, 'samtools', 'tab_dna'),
        script_p='my_samtools_SNP_pipeline',
        mindepth= '0'
    conda:
        ENV_FILES_POOL.find_yaml('old_mapping')
    threads: 1
    shell:
        '''     
        ln -fs $(pwd)/{input.STAMPY_SAM} $(pwd)/{params.prefix}.sam
        {params.script_p} {params.prefix} {input.REF} {params.mindepth}
        '''
