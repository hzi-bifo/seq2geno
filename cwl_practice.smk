rule samtools_sort:
    input:
        input="mapped/{sample}.bam"
    output:
        output_name="mapped/{sample}.sorted.bam"
    params:
        threads=lambda wildcards, threads: threads,
        memory="4G"
    threads: 8
    cwl:
        "https://github.com/common-workflow-language/workflows/blob/"
        "fb406c95/tools/samtools-sort.cwl"
