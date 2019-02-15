rule covert_nonsyn_snps_table:
    input:
        os.path.join(TMP_D, 'non-syn_SNPs_aa.tab')
    output: 
        NONSYN_SNPS_OUT
    params:
        strains=DNA_READS.index.values.tolist()
    script: 'collect_snps_data.py'

rule covert_all_snps_table:
    input:
        os.path.join(TMP_D, 'all_SNPs_aa.tab')
    output: 
        ALL_SNPS_OUT
    params:
        strains=DNA_READS.index.values.tolist()
    script: 'collect_snps_data.py'

