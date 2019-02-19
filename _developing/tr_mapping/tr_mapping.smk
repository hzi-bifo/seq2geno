rule tr_mapping:
    input:
        REF=REF_FA,
        FQ=lambda wildcards: DNA_READS[wildcards.strain],
        REF_INDEX=REF_FA+".bwt"
    output:
        bwa_bam="{TMP_D}/{strain}/bwa/tr_bwa.raw.bam",
        bwa_sorted_bam=temp("{TMP_D}/{strain}/bwa/tr_bwa.sorted.bam"),
        bwa_nonrdnt_bam=temp("{TMP_D}/{strain}/bwa/tr_bwa.bam")
    threads: 10
    shell:
        """
        bwa mem -v 2 -M -t {params.cores} {input.REF} \
{input.FQ} |\
samtools view -b -@ {threads} - \
> {output.bwa_raw_bam}
        bamtools sort -in {output.bwa_raw_bam} -out {output.bwa_sorted_bam}
        samtools rmdup -S {output.bwa_sorted_bam} {output.bwa_nonrdnt_bam}
        """
