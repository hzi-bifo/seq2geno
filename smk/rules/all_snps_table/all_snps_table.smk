rule all_snps_table:
    input:
        snps_list=temp(os.path.join(TMP_D, 'DNA_Pool1.tab')),
        ref_gbk=config['reference_annotation']
    output:
        all_output= os.path.join(TMP_D, 'all_SNPs.tab')
    params:
        script_f='lib/snps/Snp2Amino.py'
    shell:
        """
        source activate py27
        python {params.script_f} -f {input.snps_list} -g {input.ref_gbk} \
-n all \
-o {output.all_output} 
        source deactivate
        """

