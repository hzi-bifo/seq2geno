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
<<<<<<< HEAD
        strains= DNA_READS.index.values.tolist(),
=======
        strains= STRAINS,
>>>>>>> 9a517d402eeab44360883bb68c055d5209aee4c2
        bcftools_bin= 'bcftools'
    script: 'lib/vcf_to_feat_table.py'
