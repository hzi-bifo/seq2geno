rule convert_to_flatcount:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/stampy/rna.sam'
    output:
        FLTCNT=temp('{TMP_D}/{strain}/rna.flatcount')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params.SAM2ART} -f -s 2 {input.STAMPY_SAM}  > {output} '

