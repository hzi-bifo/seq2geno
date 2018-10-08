rule ng_mask_coding_regions:
    input:
        all_vcf_gz='{TMP_D}/freebayes/multisample.vcf.gz',
        coding_bed_out=temp('{TMP_D}/freebayes/ref.coding.bed')
    output:
        coding_vcf_gz="{TMP_D}/freebayes/multisample.vcf.coding.gz",
        coding_vcf_gz_index="{TMP_D}/freebayes/multisample.vcf.coding.gz.tbi"
    params:
        bedtools_bin= 'bedtools intersect',
        bgzip_bin= 'bgzip',
        tabix_bin= 'tabix'
    shell:
        """
        {params.bedtools_bin} -a {input.all_vcf_gz} -b {input.coding_bed_out}\
        -header -wa -u | {params.bgzip_bin} -c > {output.coding_vcf_gz}
        {params.tabix_bin} -p vcf {output.coding_vcf_gz}
        """

