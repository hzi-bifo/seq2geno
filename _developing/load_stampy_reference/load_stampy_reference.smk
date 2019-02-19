rule load_stampy_reference:
    input:
        REF_FA=REF_FA
    output:
        temp(REF_FA+".stidx"),
        temp(REF_FA+".sthash")
    params:
        STAMPY=STAMPY_EXE
    conda:
        ENV_FILES_POOL.find_yaml('py27')
    shell:
        """
        {params.STAMPY} -G {input} {input}
        {params.STAMPY} -g {input} -H {input}
        """

