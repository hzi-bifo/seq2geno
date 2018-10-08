rule gene_counts_paired_by_salmon:
    input:
        FQ1=lambda wildcards: RNA_READS[wildcards.strain][0],
        FQ2=lambda wildcards: RNA_READS[wildcards.strain][1],
        REF_SALMON_INDEX_DIR= "{TMP_D}/salmon_index"
    output:
        SALMON_RAW_OUTPUT= directory('{TMP_D}/{strain}/salmon.2')
    params:
        cores=CORES,
        salmon_bin= SOFTWARE['epr_quantifior']
    shell:
        """
        source activate salmon_env
        {params.salmon_bin} quant -i {input.REF_SALMON_INDEX_DIR} \
    --gcBias \
    --sigDigits 0 \
    -l A \
    -1 {input.FQ1} \
    -2 {input.FQ2} \
    -p {params.cores} \
    -o {output.SALMON_RAW_OUTPUT}
        source deactivate
        """

