rule load_reference:
    input:
        REF_FA=REF_FA
    output:
        temp(REF_FA+".fai"),
        temp(REF_FA+".bwt")
    conda:
        ENV_FILES_POOL.find_yaml('old_mapping')
    shell:
        """
        samtools faidx {input.REF_FA}
        bwa index {input}
        """

