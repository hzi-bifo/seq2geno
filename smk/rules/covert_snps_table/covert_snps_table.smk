rule covert_snps_table:
    input:
        non_syn_output= os.path.join(TMP_D, 'non-syn_SNPs.tab'),
        syn_output= os.path.join(TMP_D, 'all_SNPs.tab')
    output: 
        config['nonsyn_snps_table'],
        config['all_snps_table']
    params:
        strains= STRAINS
    script: 'collect_snps_data.py'

