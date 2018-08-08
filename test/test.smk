
rule all:
    input:
        'aln.list'

rule last:
    output:
        'aln.list'
    input:
        dynamic('{fam}.aln')
    shell:
        """
        echo \'{input}\'| sed 's/ /\\n/g' >> {output}
        """

rule middle_to_expand:
    output:
        '{fam}.aln'
    input:
        '{fam}.fa'
    shell:
        """
        cat {input} | sed \'s/gene/family/\' > {output}
        """

rule source_information:
    output:
        dynamic('{fam}.fa')
    input:
        'fam.list'
    shell:
        """
        parallel -j 1 \'echo {{}} > {{}}.fa\' \
        ::: `cat {input}`
        """
