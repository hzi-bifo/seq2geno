#####
### can try wrappers
rule create_assembly:
    input:
        READS1=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads1'],
        READS2=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads2']
    output:
        assembly="{TMP_D}/{strain}/{assembler}/scaffolds.fasta"
    params:
        cores=CORES,
#        assembler_outdir="{TMP_D}/{assembler}/{strain}",
        assembler_out="scaffolds.fasta",
        SPADES_OPT='--careful',
        SPADES_BIN='spades.py'
    shell:
        """
        {params.SPADES_BIN} {params.SPADES_OPT} \
        --threads {params.cores} \
        -o {wildcards.TMP_D}/{wildcards.strain}/{wildcards.assembler} \
        -1 {input[READS1]} -2 {input[READS2]}
        """
