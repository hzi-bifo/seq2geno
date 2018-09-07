'''
Estimate the number of reads using Salmon
'''

rule create_salmon_ref_index:
    input:
        REF_SALMON_INDEX_INPUT= temp("{TMP_D}/reference.cds.fa")
    output:
        REF_SALMON_INDEX_DIR= temp(directory("{TMP_D}/salmon_index"))
    params: 
        salmon_bin= 'salmon'
    shell:
        """
        source activate salmon_env
        {params.salmon_bin} index -t {input.REF_SALMON_INDEX_INPUT} -i \
        {output.REF_SALMON_INDEX_DIR}
        source deactivate
        """

rule for_salmon_extract_cds_seqs:
    input:
        ref_gbk=REF_GBK
    output: 
        REF_SALMON_INDEX_INPUT= temp("{TMP_D}/reference.cds.fa")
    params:
        fea_type='gene',
        name_field= 'locus_Tag'
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

rule gene_counts_by_salmon:
    input:
        FQ=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'rna_reads'],
        REF_SALMON_INDEX_DIR= temp(directory("{TMP_D}/salmon_index"))
    output:
        SALMON_OUTPUT= directory('{TMP_D}/{strain}/salmon')
    params:
        cores=CORES,
        salmon_bin= 'salmon'
    shell:
        """
        source activate salmon_env
        {params.salmon_bin} quant -i {input.REF_SALMON_INDEX_DIR} \
        --gcBias \
        -l A -r {input.FQ} \
        -p {params.cores} \
        -o {output.SALMON_OUTPUT}
        source deactivate
        """

rule create_and_make_expr_table:
    input :
        SALMON_OUTPUTS= expand(directory('{TMP_D}/{strain}/salmon'),
TMP_D=TMP_D, strain= STRAINS)
    output:
        expr_table=config['expr_table']
    params:
        strains= STRAINS,
        tmp_d= TMP_D,
        salmon_raw_output= 'quant.sf'
    run:
        import pandas as pd
        import os

        files_dict= {s: os.path.join(params.tmp_d, s, 'salmon', 'quant.sf') for
s in params.strains}
        ## open and read the Salmon outputs
        rnum_dict= {s: pd.read_csv(files_dict[s], sep= '\t', header= 0,
            index_col= 0, dtype= str)['NumReads'] for s in params.strains}
        ## convert to a matrix
        rnum_df= pd.DataFrame(data= rnum_dict).transpose() # strains in rows
        rnum_df.to_csv(output.expr_table, 
            sep= '\t', 
            na_rep= 'NA', 
            header=True, index= True)

        
