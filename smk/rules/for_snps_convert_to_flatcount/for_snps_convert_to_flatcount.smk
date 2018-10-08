rule for_snps_convert_to_flatcount:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/stampy/dna_for_tab.sam'
    output:
        FLTCNT=temp('{TMP_D}/{strain}/dna_for_tab.flatcount')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params.SAM2ART} -f -s 2 -p {input.STAMPY_SAM}  > {output} '

