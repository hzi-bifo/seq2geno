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

