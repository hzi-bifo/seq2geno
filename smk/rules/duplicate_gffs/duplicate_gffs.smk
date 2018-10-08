rule duplicate_gffs:
    ## merely for roary
    input:
        gff_output="{TMP_D}/{strain}/{assembler}/{annotator}/de_novo.gff"
    output:
        gff_output_copy=temp("{TMP_D}/{strain}/{assembler}/{annotator}/{strain}.gff")
    shell:
        """
        ln -fs {input[gff_output]} {output[gff_output_copy]}
        """
        
