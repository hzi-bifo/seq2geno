rule spades_single_create_assembly:
    input:
        READS=lambda wildcards: DNA_READS[wildcards.strain],
    output:
        spades_raw_out_dir=directory("{TMP_D}/{strain}/spades.1")
    params:
        assembler_out="scaffolds.fasta",
        SPADES_OPT='--careful',
        SPADES_BIN='spades.py'
    threads: 1
    shell:
        """
        {params.SPADES_BIN} {params.SPADES_OPT} \
--threads {threads} \
-o {output.spades_raw_out_dir} \
-s {input.READS}
        """
