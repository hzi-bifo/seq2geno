rule covert_snps_table:
    input:
        non_syn_output= os.path.join(TMP_D, 'non-syn_SNPs_aa.tab'),
        all_output= os.path.join(TMP_D, 'all_SNPs_aa.tab')
    output: 
        NONSYN_SNPS_OUT,
        SYN_SNPS_OUT
    params:
        strains=DNA_READS.index.values.tolist()
    script: 'collect_snps_data.py'

