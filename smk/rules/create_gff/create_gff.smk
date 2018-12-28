rule create_gff:
    input:
        assembly="{TMP_D}/{strain}/{assembler}/scaffolds.fasta"
    output:
        ffn_output="{TMP_D}/{strain}/{assembler}/prokka/de_novo.ffn",
        gff_output="{TMP_D}/{strain}/{assembler}/prokka/de_novo.gff"
    params:
        anno_prefix='de_novo',
        outdir= lambda wildcards: os.path.join(wildcards.TMP_D,
wildcards.strain, wildcards.assembler, 'prokka') 
    threads:1
    shell:
        """
        prokka --cpus {threads} --force --prefix {params.anno_prefix} \
--locustag {wildcards.strain} \
--outdir {params.outdir} \
{input.assembly}
        """
