
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

rule compute_gpa_raw_table:
    input:
        gffs=lambda wildcards: [os.path.join(wildcards.TMP_D, strain,
SOFTWARE['assembler'], SOFTWARE['annotator'], strain+'.gff') for strain in
DNA_READS.index.values.tolist()]
    output:
        roary_gpa="{TMP_D}/roary/gene_presence_absence.csv",
        roary_gpa_rtab='{TMP_D}/roary/gene_presence_absence.Rtab'
    params:
        ROARY_BIN='roary',
    threads:1
    conda: ENV_FILES_POOL.find_yaml('roary_env')
    shell:
        ## remove the roary folder created by snakemake first. 
        ## Otherwise, roary would create another and put all the output files in another automatically created folder
        '''
        echo $PERL5LIB
        rm -r {wildcards.TMP_D}/roary
        {params.ROARY_BIN} \
-f {wildcards.TMP_D}/roary -v -p \
{threads} -g 100000 {input.gffs}
        '''
rule:

rule:
