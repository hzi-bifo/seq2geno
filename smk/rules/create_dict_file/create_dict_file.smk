rule create_dict_file: 
    input: 
        flt_vcf_files=expand("{strain}.flt.vcf", 
            strain=STRAINS),
        flatcount_files=expand("{strain}.flatcount", 
            strain=STRAINS)
    output:
        dict_file= temp(os.path.join(TMP_D, 'dict.txt'))
    params:
        strains=STRAINS
    shell:
        """
        echo {input.flt_vcf_files}| \
        sed 's/\.flt\.vcf\W*/\\n/g'| \
        grep '\w' > {output.dict_file}
        """


