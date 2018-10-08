rule make_assembly:
    input:
        reads1
        reads2
    output:
        scaf_f=TMP_D+'/assemblies/{strain}/scaffolds.fasta'
    shell:
        "spades.py --careful --threads {CORES} -o {TMP_D}/assemblies/{wildcards.strain} -1 {input[reads1]} -2 {input[reads2]}"
'''
