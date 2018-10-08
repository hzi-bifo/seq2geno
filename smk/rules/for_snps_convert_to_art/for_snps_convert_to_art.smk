rule for_snps_convert_to_art:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/stampy/dna_for_tab.sam'
    output:
        ART=temp('{TMP_D}/{strain}/dna_for_tab.art')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params.SAM2ART} -s 2 -p -4 {input.STAMPY_SAM} > {output} '

