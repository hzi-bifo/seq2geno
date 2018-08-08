#####
### can try wrappers
rule create_assembly:
    input:
        READS1=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads1'],
        READS2=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads2']
    output:
        assembly="{TMP_D}/{strain}/{assembler}/assem.fa"
    params:
        cores=CORES,
#        assembler_outdir="{TMP_D}/{assembler}/{strain}",
        assembler_out="scaffolds.fasta"
    shell:
        """
        spades.py --careful --threads {params.cores} -o {wildcards.TMP_D}/{wildcards.strain}/{wildcards.assembler} \
        -1 {input[READS1]} -2 {input[READS2]}
        cp {wildcards.TMP_D}/{wildcards.strain}/{wildcards.assembler}/{params.assembler_out} \
        {output[assembly]}
        """
'''
rule load_reads:
    output:
        READS1=
        READS2=
'''
