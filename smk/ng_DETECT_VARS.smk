'''
Purpose:
Call population variants with freebayes

Output:
    multisample_vcf_gz="{TMP_D}/freebayes/multisample.vcf.gz"),
'''
    
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

rule ng_sort_bam:
    input:
        paired_bam="{TMP_D}/{strain}/stampy/remap.bam"
    output:
        sorted_bam=temp("{TMP_D}/{strain}/stampy/remap_sorted.bam"),
        sorted_bam_index=temp("{TMP_D}/{strain}/stampy/remap_sorted.bam.bai")
    params:
        bamtools_bin='bamtools',
        samtools_bin='samtools'
    shell:
        """
        bamtools sort -in {input} -out {output[sorted_bam]}
        samtools index {output[sorted_bam]}
        """

rule ng_sam2bam:
    input:
        STAMPY_SAM= temp('{TMP_D}/{strain}/stampy/remap.sam')
    output:
        STAMPY_BAM= temp('{TMP_D}/{strain}/stampy/remap.bam')
    params:
        CORES=CORES
    shell:
        """
        source activate Ariane_dna
        samtools view -bS -@ {params.CORES} \
        {input.STAMPY_SAM} > {output.STAMPY_BAM}
        source deactivate
        """

rule ng_stampy_remapping:
    input:
        REF_STAMPY_INDEX=REF_FA+".stidx",
        REF_INDEXDICT=REF_FA+".sthash",
        BWA_BAM= '{TMP_D}/{strain}/stampy/pre_bwa.bam'
    output:
        STAMPY_SAM= temp('{TMP_D}/{strain}/stampy/remap.sam')
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY=STAMPY_EXE
    shell:
        """
        source activate py27
        {params[STAMPY]} \
        --readgroup=ID:{wildcards.strain},SM:{wildcards.strain}\
        -g {params.REF_PREFIX} -h {params.REF_PREFIX} \
        -t{params.CORES}  --bamkeepgoodreads -M  {input.BWA_BAM}\
        > {output.STAMPY_SAM}
        source deactivate
        """

rule ng_redirect_bwa_result:
    '''
    The switch to determine which workflow (single or paired reads) to go
    '''
    input:
        BWA_RAW_BAM=
            lambda wildcards: '{}/{}/stampy/pre_bwa.{}.bam'.format(
            wildcards.TMP_D, wildcards.strain, '1' if
            (len(DNA_READS[wildcards.strain]) == 1) else '2')
    output:
        BWA_BAM= '{TMP_D}/{strain}/stampy/pre_bwa.bam'
    shell:
        '''
        mv {input.BWA_RAW_BAM} {output.BWA_BAM}
        '''

rule ng_single_read_bwa_mapping:
    input:
        FQ1=lambda wildcards: DNA_READS[wildcards.strain][0],
        REF=REF_FA,
        REF_BWA_INDEX=REF_FA+".fai",
        REF_BWA_INDEXDICT=REF_FA+".bwt"
    output:
        P_BWA_SAI1= temp('{TMP_D}/{strain}/stampy/bwa1.sai'),
        BWA_BAM= temp('{TMP_D}/{strain}/stampy/pre_bwa.1.bam')

    params:
        bwa_exe= 'bwa',
        samtools_exe= 'samtools',
        result_dir= lambda wildcards: os.path.join(wildcards.TMP_D,
wildcards.strain, 'strain'),
        BWA_OPT='-q10',
        CORES=CORES 
    shell:
        """
        bwa aln {params.BWA_OPT} -t{params.CORES} {input[REF]} {input[FQ1]} \
        > {output[P_BWA_SAI1]}

        bwa samse -r '@RG\\tID:{wildcards.strain}\\tSM:{wildcards.strain}' \
        {input.REF} {output.P_BWA_SAI1} \
        {input.FQ1} | \
        samtools view -bS -@ {params.CORES} \
        > {output.BWA_BAM}
        """

rule ng_paired_read_bwa_mapping:
    input:
        FQ1=lambda wildcards: DNA_READS[wildcards.strain][0],
        FQ2=lambda wildcards: DNA_READS[wildcards.strain][1],
        REF=REF_FA,
        REF_BWA_INDEX=REF_FA+".fai",
        REF_BWA_INDEXDICT=REF_FA+".bwt"
    output:
        P_BWA_SAI1= temp('{TMP_D}/{strain}/stampy/bwa1.sai'),
        P_BWA_SAI2= temp('{TMP_D}/{strain}/stampy/bwa2.sai'),
        BWA_BAM= temp('{TMP_D}/{strain}/stampy/pre_bwa.2.bam')

    params:
        bwa_exe= 'bwa',
        samtools_exe= 'samtools',
        result_dir= lambda wildcards: os.path.join(wildcards.TMP_D,
wildcards.strain, 'strain'),
        BWA_OPT='-q10',
        CORES=CORES 
    shell:
        """
        bwa aln {params.BWA_OPT} -t{params.CORES} {input[REF]} {input[FQ1]} \
        > {output[P_BWA_SAI1]}

        bwa aln {params.BWA_OPT} -t{params.CORES} {input[REF]} {input[FQ2]} \
        > {output[P_BWA_SAI2]} 

        bwa sampe -r '@RG\\tID:{wildcards.strain}\\tSM:{wildcards.strain}' \
        {input.REF} {output.P_BWA_SAI1} {output.P_BWA_SAI2} \
        {input.FQ1} {input.FQ2} | \
        samtools view -bS -@ {params.CORES} \
        > {output.BWA_BAM}
        """

