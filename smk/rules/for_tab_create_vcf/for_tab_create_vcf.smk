rule for_tab_create_vcf:
    input:
        REF=REF_FA,
        REF_FA_INDEX=REF_FA+".fai",
        sorted_bam= '{TMP_D}/{strain}/stampy/tab_dna.bam',
        sorted_bam_index="{TMP_D}/{strain}/stampy/tab_dna.bam.bai"
    output:
        bcf=temp('{TMP_D}/{strain}/samtools/tab_dna.raw.bcf')
    params: 
        minDepth= 0,
        CORES=CORES
    shell:
        """
        source activate old_mapping
        samtools mpileup -uf {input.REF} {input.sorted_bam} |\
bcftools view -bvcg \
- > {output.bcf}
        source deactivate
        """ 

