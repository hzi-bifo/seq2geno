'''
Purpose:
Convert the variant calling results from vcf files to (binary) features tables

Output:
        syn_snps_table= config['syn_snps_table'],
        nonsyn_snps_table= config['nonsyn_snps_table']
'''
rule vcf_to_snp_tables:
    input:
        coding_vcf_gz=os.path.join(TMP_D, 'freebayes',
'multisample.vcf.coding.gz'),
        igr_vcf_gz= os.path.join(TMP_D, 'freebayes',
'multisample.vcf.igr.gz'),
        ref_gbk=config['reference_annotation']
    output:
        syn_snps_table= config['syn_snps_table'],
        nonsyn_snps_table= config['nonsyn_snps_table']
    params:
        strains= STRAINS,
        bcftools_bin= 'bcftools'
    script: 'lib/vcf_to_feat_table.py'
