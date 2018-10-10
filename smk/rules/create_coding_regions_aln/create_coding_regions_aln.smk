rule create_coding_regions_aln:
    ## Sort sequences by gene family, align each family, and concatenate the
    ## consensus sequence alignments
    input:
        cons_coding_seqs_every_strain=expand(
            "{TMP_D}/{strains}/{caller}/cons.fa",
            TMP_D= TMP_D, 
            strains= DNA_READS.index.values.tolist(), 
            caller= 'freebayes')
    output:
        one_big_aln='{TMP_D}/OneBig.aln'

    params:
        CORES=CORES,
        TMP_D=TMP_D+"/families",
        STRAINS=DNA_READS.index.values.tolist()
    script: 'makeAlignment.py'
