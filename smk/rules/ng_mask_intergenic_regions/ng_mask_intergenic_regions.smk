rule ng_mask_intergenic_regions:
    input:
        all_vcf_gz='{TMP_D}/freebayes/multisample.vcf.gz',
        coding_bed_out=temp('{TMP_D}/freebayes/ref.coding.bed')
    output:
        igr_vcf_gz= "{TMP_D}/freebayes/multisample.vcf.igr.gz",
        igr_vcf_gz_index= "{TMP_D}/freebayes/multisample.vcf.igr.gz.tbi"
    params:
        bedtools_bin= 'bedtools intersect',
        bgzip_bin= 'bgzip',
        tabix_bin= 'tabix'
    shell:
        """
        {params.bedtools_bin} -a {input.all_vcf_gz} -b {input.coding_bed_out}\
        -header -wa -v | {params.bgzip_bin} -c > {output.igr_vcf_gz}
        {params.tabix_bin} -p vcf {output.igr_vcf_gz}
        """

