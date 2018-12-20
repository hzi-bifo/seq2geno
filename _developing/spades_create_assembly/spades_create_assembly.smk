rule spades_paired_create_assembly:
    input:
        READS=lambda wildcards: DNA_READS[wildcards.strain]
    output:
        spades_raw_out_dir=directory("{TMP_D}/{strain}/spades.2")
    params:
        assembler_out="scaffolds.fasta",
        SPADES_OPT='--careful',
        SPADES_BIN='spades.py'
    threads: 1
    conda: 
    shell:
        """
        n=$(echo {input.READS}| awk -F' ' '{{print NF}}') 
        echo $n
        if test $n -eq 2  
        then
            {params.SPADES_BIN} {params.SPADES_OPT} \
--threads {threads} \
-o {output.spades_raw_out_dir} \
-1 {input[READS1]} -2 {input[READS2]}
        """

