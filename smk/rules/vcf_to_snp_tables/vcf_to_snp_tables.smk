rule vcf_to_snp_tables:
    input:
        coding_vcf_gz=os.path.join(TMP_D, 'freebayes',
'multisample.vcf.coding.gz'),
        igr_vcf_gz= os.path.join(TMP_D, 'freebayes',
'multisample.vcf.igr.gz'),
        ref_gbk=REF_GBK
    output:
        syn_snps_table= SYN_SNPS_OUT,
        nonsyn_snps_table= NONSYN_SNPS_OUT
    params:
        strains= DNA_READS.index.values.tolist(),
        bcftools_bin= 'bcftools'
    script: 'vcf_to_feat_table.py'
