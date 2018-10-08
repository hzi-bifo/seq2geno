rule for_salmon_create_ref_index:
    input:
        REF_SALMON_INDEX_INPUT= "{TMP_D}/reference.cds.fa"
    output:
        REF_SALMON_INDEX_DIR= temp(directory("{TMP_D}/salmon_index"))
    params: 
        salmon_bin= SOFTWARE['epr_quantifior']
    shell:
        """
        source activate salmon_env
        {params.salmon_bin} index -t {input.REF_SALMON_INDEX_INPUT} -i \
        {output.REF_SALMON_INDEX_DIR}
        source deactivate
        """

