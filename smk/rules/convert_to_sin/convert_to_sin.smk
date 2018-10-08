rule convert_to_sin:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/stampy/rna.sam'
    output:
        SIN=temp('{TMP_D}/{strain}/rna.sin')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params.SAM2ART} -s 2 -l {input.STAMPY_SAM}  > {output} '

