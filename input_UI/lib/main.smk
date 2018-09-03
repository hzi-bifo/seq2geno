rule all:
    input:
        'tmp.txt'

rule some_process:
    input:
        'files.txt'
    output:
        'tmp.txt'
    shell:
        """
        cat {input} | head -n 5 > {output}
        """
