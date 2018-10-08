rule tr_single_mapping:
    input:
        REF=REF_FA,
        FQ1=lambda wildcards: DNA_READS[wildcards.strain][0],
        REF_INDEX=REF_FA+".bwt"
    output:
        bwa_raw_bam=temp("{TMP_D}/{strain}/bwa/tr_bwa.raw.1.bam"),
        bwa_sorted_bam=temp("{TMP_D}/{strain}/bwa/tr_bwa.sorted.1.bam"),
        bwa_nonrdnt_bam=temp("{TMP_D}/{strain}/bwa/tr_bwa.1.bam")
    params:
        cores= CORES
    shell:
        """
        bwa mem -v 2 -M -t {params.cores} {input.REF} \
{input.FQ1} |\
samtools view -b -@ {params.cores} - \
> {output.bwa_raw_bam}
        bamtools sort -in {output.bwa_raw_bam} -out {output.bwa_sorted_bam}
        samtools rmdup -s {output.bwa_sorted_bam} {output.bwa_nonrdnt_bam}
        """

