rule create_gff:
    input:
        assembly="{TMP_D}/{strain}/{assembler}/scaffolds.fasta"
    output:
        ffn_output="{TMP_D}/{strain}/{assembler}/{annotator}/de_novo.ffn",
        gff_output="{TMP_D}/{strain}/{assembler}/{annotator}/de_novo.gff"
    params:
        anno_prefix='de_novo',
        cores=CORES
    shell:
        """
        prokka --cpus {params.cores} --force --prefix {params.anno_prefix} \
--locustag {wildcards.strain} \
--outdir \
{wildcards.TMP_D}/{wildcards.strain}/{wildcards.assembler}/{wildcards.annotator} \
{input.assembly}
        """
