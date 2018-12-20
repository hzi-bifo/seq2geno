rule all_snps_table:
    input:
        snps_list=temp(os.path.join(TMP_D, 'all_SNPs.tab')),
        ref_gbk=REF_GBK
    output:
        all_output= os.path.join(TMP_D, 'all_SNPs_aa.tab')
    params:
#        script_f='Snp2Amino.edit.py'
        script_f='Snp2Amino.py'
    conda: ENV_FILES_POOL.find_yaml('old_mapping')
    threads: 1
    shell:
        """
        {params.script_f} -f {input.snps_list} -g {input.ref_gbk} \
-n all \
-o {output.all_output} 
        """

