rule rna_redirect_stampy_result:
    input:
        STAMPY_RAW_SAM=
            lambda wildcards: '{}/{}/stampy/rna.{}.sam'.format(
            wildcards.TMP_D, wildcards.strain, '1' if
            (len(RNA_READS[wildcards.strain]) == 1) else '2')
    output:
        STAMPY_SAM='{TMP_D}/{strain}/stampy/rna.sam'
    shell:
        '''
        mv {input.STAMPY_RAW_SAM} {output.STAMPY_SAM}
        '''

