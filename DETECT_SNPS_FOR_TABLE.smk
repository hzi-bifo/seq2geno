#####
# Two-stage mapping with stampy
'''
rule for_tab_index_vcf:
    input:
        vcf_gz="{TMP_D}/{strain}/{mapper}/st_vcf.gz"
    output:
        vcf_gz_index= "{TMP_D}/{strain}/{mapper}/st_vcf.gz.tbi"
    shell:
        """
        tabix -p vcf {input[vcf_gz]}
        """
'''
rule for_tab_create_vcf:
    input:
        REF=REF_FA,
        REF_FA_INDEX=REF_FA+".fai",
        BAM="{TMP_D}/{strain}/{mapper}/st_sorted.bam",
        BAM_INDEX="{TMP_D}/{strain}/{mapper}/st_sorted.bam.bai"
    output:
        bcf_out=temp("{TMP_D}/{strain}/{mapper}/st_vcf.bcf")
        vcf_out=temp("{TMP_D}/{strain}/{mapper}/st_vcf.flt.vcf")
#        vcf_gz="{TMP_D}/{strain}/{mapper}/st_vcf.gz"
    params: 
        VCFUTIL_EXE='lib/vcfutils.pl',
        minDepth= 0,
        CORES=CORES
        
    shell:
        """
        samtools mpileup -uf {input.REF} {input.BAM} | bcftools view -bvcg - > {output.bcf_out}
        bcftools view {output.bcf_out} | {params.VCFUTILS} varFilter -d {params.minDepth} > {output.vcf_out}
        """ 

rule for_tab_sort_bam:
    input:
        paired_bam="{TMP_D}/{strain}/{mapper}/st_paired.bam"
    output:
        sorted_bam=temp("{TMP_D}/{strain}/{mapper}/st_sorted.bam"),
        sorted_bam_index=temp("{TMP_D}/{strain}/{mapper}/st_sorted.bam.bai")
    shell:
        """
        bamtools sort -in {input} -out {output[sorted_bam]}
        samtools index {output[sorted_bam]}
        """

rule for_tab_sam2bam:
    input:
        STAMPY_SAM= temp('{TMP_D}/{strain}/{mapper}/st_paired.sam')
    output:
        STAMPY_BAM= temp('{TMP_D}/{strain}/{mapper}/st_paired.bam')
    params:
        CORES=CORES
    shell:
        """
        samtools view -bS -@ {params.CORES} \
        {input.STAMPY_SAM} > {output.STAMPY_BAM}
        """

rule for_tab_stampy_remapping:
    input:
        REF_STAMPY_INDEX=REF_FA+".stidx",
        REF_INDEXDICT=REF_FA+".sthash",
        BWA_BAM= '{TMP_D}/{strain}/{mapper}/dna.bwa.bam'
    output:
        STAMPY_SAM= temp('{TMP_D}/{strain}/{mapper}/st_paired.sam')
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY=STAMPY_EXE
    shell:
        """
        source activate py27
        {params[STAMPY]} \
        -g {params.REF_PREFIX} -h {params.REF_PREFIX} \
        -t{params.CORES}  --bamkeepgoodreads -M  {input[BWA_BAM]}\
        > {output.STAMPY_SAM}
        """

rule for_tab_paired_read_bwa_mapping:
    input:
        FQ1=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads1'],
        FQ2=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads2'],
        REF=REF_FA,
        REF_BWA_INDEX=REF_FA+".fai",
        REF_BWA_INDEXDICT=REF_FA+".bwt"
    output:
        P_BWA_SAI1= temp('{TMP_D}/{strain}/{mapper}/dna1.bwa.sai'),
        P_BWA_SAI2= temp('{TMP_D}/{strain}/{mapper}/dna2.bwa.sai'),
        BWA_BAM= temp('{TMP_D}/{strain}/{mapper}/dna.bwa.bam')

    params:
        BWA_OPT='-q10',
        CORES=CORES 
    shell:
        """
        bwa aln {params.BWA_OPT} -t{params.CORES} {input[REF]} {input[FQ1]} \
        > {output[P_BWA_SAI1]}

        bwa aln {params.BWA_OPT} -t{params.CORES} {input[REF]} {input[FQ2]} \
        > {output[P_BWA_SAI2]} 

        bwa sampe {input[REF]} {output[P_BWA_SAI1]} {output[P_BWA_SAI2]} \
        {input[FQ1]} {input[FQ2]} | \
        samtools view -bS -@ {params.CORES} \
        > {output[BWA_BAM]}
        """
'''
rule for_tab_load_reference:
    input:
        REF_FA=REF_FA
    output:
        temp(REF_FA+".fai"),
        temp(REF_FA+".bwt"),
        temp(REF_FA+".stidx"),
        temp(REF_FA+".sthash")
    params:
        STAMPY=STAMPY_EXE
    shell:
        """
        samtools faidx {input[REF_FA]}
        bwa index {input}
        source activate py27
        {params[STAMPY]} -G {input} {input}
        {params[STAMPY]} -g {input} -H {input}
        source deactivate
        """
'''
