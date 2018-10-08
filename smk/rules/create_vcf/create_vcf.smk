rule create_vcf:
    input:
        REF=REF_FA,
        REF_FA_INDEX=REF_FA+".fai",
        BAM="{TMP_D}/{strain}/bwa/tr_bwa.bam",
        BAM_INDEX="{TMP_D}/{strain}/bwa/tr_bwa.bam.bai"
    output:
        vcf_gz="{TMP_D}/{strain}/freebayes/vcf.gz"
    params: 
        CORES=CORES
    shell:
        """
        freebayes-parallel <(fasta_generate_regions.py \
{input[REF_FA_INDEX]} 100000) {params[CORES]} \
-p 1 -f {input[REF]} {input[BAM]} | \
bgzip -c > {output[vcf_gz]}
        """ 

