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

