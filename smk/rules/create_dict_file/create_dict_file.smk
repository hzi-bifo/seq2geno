rule create_dict_file: 
    input: 
        flt_vcf_files=expand("{strain}.flt.vcf", 
            strain=DNA_READS.index.values.tolist()),
        flatcount_files=expand("{strain}.flatcount", 
            strain=DNA_READS.index.values.tolist())
    output:
        dict_file= temp(os.path.join(TMP_D, 'dict.txt'))
    params:
        strains=DNA_READS.index.values.tolist()
    shell:
        """
        echo {input.flt_vcf_files}| \
        sed 's/\.flt\.vcf\W*/\\n/g'| \
        grep '\w' > {output.dict_file}
        """


