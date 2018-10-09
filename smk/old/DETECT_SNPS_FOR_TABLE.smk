'''
Purpose:
Detect variants for the snps table.

Output:
        vcf=temp("{TMP_D,"^[^/]+$"}/{strain}/samtools/tab_dna.flt.vcf")
'''
    
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
        source activate Ariane_dna
        samtools mpileup -uf {input.REF} {input.sorted_bam} |\
bcftools view -bvcg \
- > {output.bcf}
        source deactivate
        """ 

rule for_tab_sam2bam:
    input:
        STAMPY_SAM= temp('{TMP_D}/{strain}/stampy/dna_for_tab.sam')
    output:
        bam= temp('{TMP_D}/{strain}/stampy/tab_dna.bam'),
        bam_index=temp("{TMP_D}/{strain}/stampy/tab_dna.bam.bai")
    params: 
        output_prefix= lambda wildcards: os.path.join(wildcards.TMP_D,
wildcards.strain, 'stampy', 'tab_dna')
    shell:
        """
        source activate Ariane_dna
        samtools view -bS \
{input.STAMPY_SAM} > {output.bam}
        samtools sort {output.bam} \
{params.output_prefix}
        samtools index {output.bam}
        source deactivate
        """

rule for_snps_convert_to_art:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/stampy/dna_for_tab.sam'
    output:
        ART=temp('{TMP_D}/{strain}/dna_for_tab.art')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params.SAM2ART} -s 2 -p -4 {input.STAMPY_SAM} > {output} '

rule for_snps_convert_to_sin:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/stampy/dna_for_tab.sam'
    output:
        SIN=temp('{TMP_D}/{strain}/dna_for_tab.sin')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params.SAM2ART} -s 2 -l -p -4 {input.STAMPY_SAM}  > {output} '

rule for_snps_convert_to_flatcount:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/stampy/dna_for_tab.sam'
    output:
        FLTCNT=temp('{TMP_D}/{strain}/dna_for_tab.flatcount')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params.SAM2ART} -f -s 2 -p {input.STAMPY_SAM}  > {output} '

rule for_tab_redirect_stampy_result:
    input:
        STAMPY_RAW_SAM=
            lambda wildcards: '{}/{}/stampy/dna_for_tab.{}.sam'.format(
            wildcards.TMP_D, wildcards.strain, '1' if
            (len(DNA_READS[wildcards.strain]) == 1) else '2')
    output:
        STAMPY_SAM='{TMP_D}/{strain}/stampy/dna_for_tab.sam'
    shell:
        '''
        mv {input.STAMPY_RAW_SAM} {output.STAMPY_SAM}
        '''        

rule for_tab_stampy_single_mapping:
    input:
        FQ1=lambda wildcards: DNA_READS[wildcards.strain][0],
        REF=REF_FA,
        REF_BWA_INDEX=REF_FA+".fai",
        REF_BWA_INDEXDICT=REF_FA+".bwt"
    output:
        STAMPY_RAW_SAM= temp('{TMP_D}/{strain}/stampy/dna_for_tab.1.sam')
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY=STAMPY_EXE
    shell:
        """
        source activate Ariane_dna
        {params.STAMPY} --bwaoptions=\"-q10 {input.REF}\" \
-g {params.REF_PREFIX} \
-h {params.REF_PREFIX} \
-M {input.FQ1} > {output.STAMPY_RAW_SAM} 
        source deactivate
        """

rule for_tab_stampy_paired_mapping:
    input:
        FQ1=lambda wildcards: DNA_READS[wildcards.strain][0],
        FQ2=lambda wildcards: DNA_READS[wildcards.strain][1],
        REF=REF_FA,
        REF_BWA_INDEX=REF_FA+".fai",
        REF_BWA_INDEXDICT=REF_FA+".bwt"
    output:
        STAMPY_RAW_SAM= temp('{TMP_D}/{strain}/stampy/dna_for_tab.2.sam')
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY=STAMPY_EXE
    shell:
        """
        source activate Ariane_dna
        {params.STAMPY} --bwaoptions=\"-q10 {input.REF}\" \
-g {params.REF_PREFIX} \
-h {params.REF_PREFIX} \
-M {input.FQ1} {input.FQ2} > {output.STAMPY_RAW_SAM} 
        source deactivate
        """
