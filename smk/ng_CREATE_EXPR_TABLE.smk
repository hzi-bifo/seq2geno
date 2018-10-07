'''
Purpose:
Estimate the number of reads using Salmon

Output:
        expr_table=config['expr_table']
'''

rule create_and_make_expr_table:
    input :
        SALMON_OUTPUTS= expand(directory('{TMP_D}/{strain}/salmon'),
TMP_D=TMP_D, strain= RNA_READS.index.values.tolist())
    output:
        expr_table=config['expr_table']
    params:
        strains= RNA_READS.index.values.tolist(),
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

rule gene_counts_paired_by_salmon:
    input:
        FQ1=lambda wildcards: RNA_READS[wildcards.strain][0],
        FQ2=lambda wildcards: RNA_READS[wildcards.strain][1],
        REF_SALMON_INDEX_DIR= "{TMP_D}/salmon_index"
    output:
        SALMON_RAW_OUTPUT= directory('{TMP_D}/{strain}/salmon.2')
    params:
        cores=CORES,
        salmon_bin= SOFTWARE['epr_quantifior']
    shell:
        """
        source activate salmon_env
        {params.salmon_bin} quant -i {input.REF_SALMON_INDEX_DIR} \
    --gcBias \
    --sigDigits 0 \
    -l A \
    -1 {input.FQ1} \
    -2 {input.FQ2} \
    -p {params.cores} \
    -o {output.SALMON_RAW_OUTPUT}
        source deactivate
        """

rule gene_counts_single_by_salmon:
    input:
        FQ=lambda wildcards: RNA_READS[wildcards.strain][0],
        REF_SALMON_INDEX_DIR= "{TMP_D}/salmon_index"
    output:
        SALMON_RAW_OUTPUT= directory('{TMP_D}/{strain}/salmon.1')
    params:
        cores=CORES,
        salmon_bin= SOFTWARE['epr_quantifior']
    shell:
        """
        source activate salmon_env
        {params.salmon_bin} quant -i {input.REF_SALMON_INDEX_DIR} \
    --gcBias \
    -l A -r {input.FQ} \
    --sigDigits 0 \
    -p {params.cores} \
    -o {output.SALMON_RAW_OUTPUT}
        source deactivate
        """

rule for_salmon_create_ref_index:
    input:
        REF_SALMON_INDEX_INPUT= "{TMP_D}/reference.cds.fa"
    output:
        REF_SALMON_INDEX_DIR= temp(directory("{TMP_D}/salmon_index"))
    params: 
        salmon_bin= SOFTWARE['epr_quantifior']
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

        
