rule for_salmon_extract_cds_seqs:
    input:
        ref_gbk=REF_GBK
    output: 
        REF_SALMON_INDEX_INPUT= temp("{TMP_D}/reference.cds.fa")
    params:
        fea_type='gene',
        name_field= 'locus_tag'
    run:
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC
        import re
        # read the gbk
        gbk_f= input.ref_gbk
        rec= SeqIO.read(open(gbk_f, 'r', encoding='windows-1252'), 'gb')
        chromosome= rec.id
        features= [fea for fea in rec.features if ((fea.type == params.fea_type)
            & ('locus_tag' in fea.qualifiers))]

        seq_records=[]
        for fea in features:
            fasta_header=fea.qualifiers['locus_tag'][0]
            start= int(re.sub('\W','',str(fea.location.start)))
            end= int(re.sub('\W','', str(fea.location.end)))

            seq= rec.seq[start:end]
            seq_rec= SeqRecord(seq, id= fasta_header, name= '', description= '')
            seq_records.append(seq_rec)
        with open(output[0], 'w') as out_fh:
            SeqIO.write(seq_records, out_fh, 'fasta')

        
