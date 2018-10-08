rule all:
    input:
        config['tree'],
        config['expr_table'],
        config['syn_snps_table'],
        config['syn_snps_table']+'_GROUPS',
        config['syn_snps_table']+'_NON-RDNT',
        config['nonsyn_snps_table'],
        config['nonsyn_snps_table']+'_GROUPS',
        config['nonsyn_snps_table']+'_NON-RDNT'

