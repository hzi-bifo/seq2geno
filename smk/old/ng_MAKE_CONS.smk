rule create_consensus_sequence:
    ## Create consensus sequences of the coding regions 
    input:
        ref_target_seqs="{TMP_D}/reference.target_regions.fa",
        coding_vcf_gz="{TMP_D}/freebayes/multisample.vcf.coding.gz",
        coding_vcf_gz_index="{TMP_D}/freebayes/multisample.vcf.coding.gz.tbi"
    output:
        cons_coding_seqs="{TMP_D}/{strain}/cons.fa"
    shell:
        """
        bcftools consensus --sample {wildcards.strain} \
-f {input.ref_target_seqs} \
        {input.coding_vcf_gz} > {output.cons_coding_seqs}
        """


rule for_consensus_create_ref_regional_seq:
    ## Which feature to parse?
    input:
        ref_gbk=REF_GBK
    output: 
        ref_target_seqs=temp("{TMP_D}/reference.target_regions.fa")
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
            start= int(re.sub('\W','',str(fea.location.start)))
            end= int(re.sub('\W','', str(fea.location.end)))
            fasta_header='{}:{}-{}'.format(chromosome, str(start+1), str(end)) 
            seq= rec.seq[start:end]
            seq_rec= SeqRecord(seq, id= fasta_header, name= '', description= '')
            seq_records.append(seq_rec)
        with open(output[0], 'w') as out_fh:
            SeqIO.write(seq_records, out_fh, 'fasta')

