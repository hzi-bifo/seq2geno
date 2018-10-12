rule spades_paired_create_assembly:
    input:
        READS1=lambda wildcards: DNA_READS[wildcards.strain],
        READS2=lambda wildcards: DNA_READS[wildcards.strain]
    output:
        spades_raw_out_dir=directory("{TMP_D}/{strain}/spades.2")
    params:
        cores=CORES,
        assembler_out="scaffolds.fasta",
        SPADES_OPT='--careful',
        SPADES_BIN='spades.py'
    shell:
        """
        {params.SPADES_BIN} {params.SPADES_OPT} \
--threads {params.cores} \
-o {output.spades_raw_out_dir} \
-1 {input[READS1]} -2 {input[READS2]}
        """