rule redirect_salmon_result:
    input:
        SALMON_RAW_OUTPUT= lambda wildcards:
            '{}/{}/salmon.{}'.format(
            wildcards.TMP_D, wildcards.strain, '1' if
            (len(RNA_READS[wildcards.strain]) == 1) else '2')
    output: 
        SALMON_OUTPUT= directory('{TMP_D}/{strain}/salmon')
    shell:
        '''
        mv {input.SALMON_RAW_OUTPUT} {output.SALMON_OUTPUT}
        '''

