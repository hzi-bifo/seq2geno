rule re_redirect_bwa_result:
    input:
        bwa_raw_bam=
            lambda wildcards: "{}/{}/bwa/tr_bwa.{}.bam".format(
            wildcards.TMP_D, wildcards.strain, '1' if
            (len(DNA_READS[wildcards.strain]) == 1) else '2')
    output:
        bwa_bam="{TMP_D}/{strain}/bwa/tr_bwa.bam"
    shell:
        '''
        mv {input.bwa_raw_bam} {output.bwa_bam}
        '''

