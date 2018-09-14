'''
Create multi-sample vcf
'''
rule ng_multi_sample_vcf:
    input:
        vcf_files=lambda wildcards: [os.path.join(wildcards.TMP_D,
strain,'stampy', 'vcf.gz') for strain in STRAINS]
    output:
        all_vcf_gz='{TMP_D}/multi_sample_vcf.gz',
        all_vcf_gz_index='{TMP_D}/multi_sample_vcf.gz.tbi'
    params:
        cores= CORES,
        bgzip_bin= 'bgzip',
        tabix_bin= 'tabix', 
        bcftools_bin= 'bcftools merge'
    shell:
        """
        {params.bcftools_bin} -m none {input.vcf_files} \
        --threads {params.cores}|\
        {params.bgzip_bin} -c > {output.all_vcf_gz}
        {params.tabix_bin} -p vcf {output.all_vcf_gz}
        """
