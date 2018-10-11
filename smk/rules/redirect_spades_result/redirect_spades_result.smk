rule redirect_spades_result:
    input:
        spades_raw_out_dir= 
            lambda wildcards: '{}/{}/spades.{}'.format(
            wildcards.TMP_D, wildcards.strain, '1' if
            (len(DNA_READS[wildcards.strain]) == 1) else '2')
    output:
        spades_out_scf= '{TMP_D}/{strain}/spades/scaffolds.fasta'
    params:
        spades_out_dir= '{TMP_D}/{strain}/spades'
        
    shell:
        '''
        mv {input.spades_raw_out_dir} {params.spades_out_dir}
        '''
