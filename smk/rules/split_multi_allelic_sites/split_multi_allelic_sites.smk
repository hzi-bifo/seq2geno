rule split_multi_allelic_sites:
    input:
        multisample_raw_vcf_gz="{TMP_D}/freebayes/multisample_raw.vcf.gz",
        multisample_raw_vcf_gz_index="{TMP_D}/freebayes/multisample_raw.vcf.gz.tbi"
    output:
        all_vcf_gz="{TMP_D}/freebayes/multisample.vcf.gz",
        all_vcf_gz_index="{TMP_D}/freebayes/multisample.vcf.gz.tbi"
    params: 
        CORES=CORES,
        tabix_bin= 'tabix',
        bcftools_bin= 'bcftools norm'
    shell:
        """
        {params.bcftools_bin} -m -both {input.multisample_raw_vcf_gz} \
        -O z -o {output.all_vcf_gz}
        {params.tabix_bin} -p vcf {output.all_vcf_gz}
        """
         
