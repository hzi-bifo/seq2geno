rule create_coding_regions_aln:
    ## Sort sequences by gene family, align each family, and concatenate the
    ## consensus sequence alignments
    input:
        cons_coding_seqs_every_strain=expand(
            "{TMP_D}/{strain}/freebayes/cons.fa",
            TMP_D= TMP_D, 
            strain= DNA_READS.index.values.tolist())
    output:
        one_big_aln='{TMP_D}/OneBig.aln'
    threads: int(CORES)
    params:
        TMP_D=TMP_D+"/families",
        STRAINS=DNA_READS.index.values.tolist()
    script: 'makeAlignment.py'
