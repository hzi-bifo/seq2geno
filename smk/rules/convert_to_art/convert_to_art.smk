rule convert_to_art:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/stampy/rna.sam'
    output:
        ART=temp('{TMP_D}/{strain}/rna.art')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params.SAM2ART} -s 2 {input.STAMPY_SAM} > {output} '

