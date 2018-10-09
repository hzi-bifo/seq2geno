rule nonsyn_snps_table:
    input:
        snps_list=temp(os.path.join(TMP_D, 'DNA_Pool1.tab')),
        ref_gbk=REF_GBK
    output: 
        non_syn_output= os.path.join(TMP_D, 'non-syn_SNPs.tab')
    params:
        script_f='lib/snps/Snp2Amino.py'
    shell:
        """
        source activate py27
        python {params.script_f} -f {input.snps_list} -g {input.ref_gbk} \
-n non-syn \
-o {output.non_syn_output} 
        source deactivate
        """

