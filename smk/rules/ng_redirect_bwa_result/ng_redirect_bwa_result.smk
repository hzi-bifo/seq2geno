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

