rule ng_detect_var_for_pop:
    input:    
        REF=REF_FA,
        REF_FA_INDEX=REF_FA+".fai",
        BAM=lambda wildcards: [os.path.join(wildcards.TMP_D, strain, 'stampy',
'remap_sorted.bam') for strain in DNA_READS.index.values.tolist()],
        BAM_INDEX= lambda wildcards: [os.path.join(wildcards.TMP_D, strain, 'stampy', 'remap_sorted.bam.bai') for strain in DNA_READS.index.values.tolist()]
    output:
        multisample_raw_vcf_gz="{TMP_D}/freebayes/multisample_raw.vcf.gz",
        multisample_raw_vcf_gz_index="{TMP_D}/freebayes/multisample_raw.vcf.gz.tbi"
    params: 
        CORES=CORES,
        bgzip_bin= 'bgzip',
        tabix_bin= 'tabix',
        frbayes_bin= 'freebayes-parallel',
        regions_generator_bin='fasta_generate_regions.py' 
    shell:
        """
        {params.frbayes_bin} <({params.regions_generator_bin} \
        {input.REF_FA_INDEX} 100000) {params.CORES} \
        -p 1 -f {input.REF} {input.BAM} | \
        {params.bgzip_bin} -c > {output.multisample_raw_vcf_gz}
        {params.tabix_bin} -p vcf {output.multisample_raw_vcf_gz}
        """

