rule for_tab_filter_vcf:
    input:
        bcf='{TMP_D}/{strain}/samtools/tab_dna.raw.bcf'
    output:
        vcf='{TMP_D}/{strain}/samtools/tab_dna.flt.vcf'
    params: 
        VCFUTIL_EXE='lib/vcfutils.pl',
        minDepth= 0
    shell:
        """
        source activate Ariane_dna
        bcftools view {input.bcf} |\
        {params.VCFUTIL_EXE} varFilter -d {params.minDepth} > {output.vcf}
        """ 
    
