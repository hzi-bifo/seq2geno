rule ng_create_coding_region_bed:
    input:
        ref_gbk=REF_GBK
    output:
        coding_bed_out=temp('{TMP_D}/freebayes/ref.coding.bed')
    run:
        from Bio import SeqIO
        rec= SeqIO.read(open(input.ref_gbk, 'r', encoding='windows-1252'), 'gb')
        chromosome= rec.id
        # filter target regions types (gene, CDS...etc)
        features= [fea for fea in rec.features if fea.type == 'CDS']
        with open(output.coding_bed_out, 'w') as coding_bed_out_h:
            rows=[create_bed3_string(fea, chromosome) for fea in features]
            rows.sort(key= lambda x:int((x.split('\t'))[1]))
            coding_bed_out_h.write('\n'.join(rows)+'\n')

