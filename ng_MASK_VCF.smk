'''
Purpose:
For the multi-sample vcf, ,ask variants by coding regions

Output:
        igr_vcf_gz= "{TMP_D}/freebayes/multisample.vcf.igr.gz",
        coding_vcf_gz="{TMP_D}/freebayes/multisample.vcf.coding.gz",
'''

def create_bed3_string(fea, chromosome):
    '''
    Purpose:
    Create the bed-style coordinate string

    Input:
    The feature object of genbank, which should include the qualifiers "name"
    and "locus_tag"
    '''
    gene_name=fea.qualifiers['name'][0] if 'name' in fea.qualifiers else '.'
    gene_id=fea.qualifiers['locus_tag'][0] if 'locus_tag' in fea.qualifiers else '.'
    ## In bed format, it's 0-based starting for starting site while 1-based for
    ## ending site, which is exactly the parsing result of biopython
    start= re.sub('\W','',str(fea.location.start))
    end= re.sub('\W','', str(fea.location.end))
    coord_str= '\t'.join([chromosome, start, end])
    return(coord_str)

rule ng_mask_intergenic_regions:
    input:
        all_vcf_gz='{TMP_D}/freebayes/multisample.vcf.gz',
        coding_bed_out=temp('{TMP_D}/freebayes/ref.coding.bed')
    output:
        igr_vcf_gz= "{TMP_D}/freebayes/multisample.vcf.igr.gz",
        igr_vcf_gz_index= "{TMP_D}/freebayes/multisample.vcf.igr.gz.tbi"
    params:
        bedtools_bin= 'bedtools intersect',
        bgzip_bin= 'bgzip',
        tabix_bin= 'tabix'
    shell:
        """
        {params.bedtools_bin} -a {input.all_vcf_gz} -b {input.coding_bed_out}\
        -header -wa -v | {params.bgzip_bin} -c > {output.igr_vcf_gz}
        {params.tabix_bin} -p vcf {output.igr_vcf_gz}
        """

rule ng_mask_coding_regions:
    input:
        all_vcf_gz='{TMP_D}/freebayes/multisample.vcf.gz',
        coding_bed_out=temp('{TMP_D}/freebayes/ref.coding.bed')
    output:
        coding_vcf_gz="{TMP_D}/freebayes/multisample.vcf.coding.gz",
        coding_vcf_gz_index="{TMP_D}/freebayes/multisample.vcf.coding.gz.tbi",
    params:
        bedtools_bin= 'bedtools intersect',
        bgzip_bin= 'bgzip',
        tabix_bin= 'tabix'
    shell:
        """
        {params.bedtools_bin} -a {input.all_vcf_gz} -b {input.coding_bed_out}\
        -header -wa -u | {params.bgzip_bin} -c > {output.coding_vcf_gz}
        {params.tabix_bin} -p vcf {output.coding_vcf_gz}
        """

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

