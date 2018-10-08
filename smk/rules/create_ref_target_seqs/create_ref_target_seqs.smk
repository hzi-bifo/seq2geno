rule create_ref_target_seqs:
    ## Which feature to parse?
    input:
        ref_gbk=REF_GBK
    output: 
        ref_target_seqs="{TMP_D}/reference.target_regions.fa"
    params:
        fea_type='CDS'
    run:
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC
        import re
        # read the gbk
        gbk_f= input[0]
        rec= SeqIO.read(open(gbk_f, 'r', encoding='windows-1252'), 'gb')
        chromosome= rec.id
        # filter target regions types (gene, CDS...etc)
        features= [fea for fea in rec.features if fea.type == params['fea_type']]

        seq_records=[]
        for fea in features:
            #print(fea.strand)
            start= int(re.sub('\W','',str(fea.location.start)))
            end= int(re.sub('\W','', str(fea.location.end)))

            seq= rec.seq[start:end]
            fasta_header='{}:{}-{}'.format(chromosome, str(start+1), str(end)) # notice that the coordinate should be gff style as the bcftools uses it
            seq_rec= SeqRecord(seq, id= fasta_header, name= '', description= '')
            seq_records.append(seq_rec)
        with open(output[0], 'w') as out_fh:
            SeqIO.write(seq_records, out_fh, 'fasta')


