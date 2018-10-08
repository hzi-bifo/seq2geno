rule index_vcf:
    input:
        vcf_gz="{TMP_D}/{strain}/freebayes/vcf.gz"
    output:
        vcf_gz_index= "{TMP_D}/{strain}/freebayes/vcf.gz.tbi"
    shell:
        """
        tabix -p vcf {input[vcf_gz]}
        """

