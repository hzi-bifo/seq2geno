rule create_dict_file: 
    input: 
        flt_vcf_files=temp(expand("{strain}.flt.vcf", 
            strain=DNA_READS.index.values.tolist())),
        flatcount_files=temp(expand("{strain}.flatcount", 
            strain=DNA_READS.index.values.tolist()))
    output:
        dict_file= temp(os.path.join(TMP_D, 'dict.txt'))
    params:
        strains=DNA_READS.index.values.tolist()
    shell:
        """
        parallel  \
 'echo {{}} \
 >> {output.dict_file}'\
 ::: {params.strains}
        """

