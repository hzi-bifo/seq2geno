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
        BAM=lambda wildcards: [os.path.join(wildcards.TMP_D, strain, 'stampy', 'remap_sorted.bam') for strain in STRAINS],
        BAM_INDEX= lambda wildcards: [os.path.join(wildcards.TMP_D, strain, 'stampy', 'remap_sorted.bam.bai') for strain in STRAINS]
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

rule ng_paired_read_bwa_mapping:
    input:
        FQ1=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads1'],
        FQ2=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads2'],
        REF=REF_FA,
        REF_BWA_INDEX=REF_FA+".fai",
        REF_BWA_INDEXDICT=REF_FA+".bwt"
    output:
        P_BWA_SAI1= temp('{TMP_D}/{strain}/stampy/bwa1.sai'),
        P_BWA_SAI2= temp('{TMP_D}/{strain}/stampy/bwa2.sai'),
        BWA_BAM= temp('{TMP_D}/{strain}/stampy/pre_bwa.bam')

    params:
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
