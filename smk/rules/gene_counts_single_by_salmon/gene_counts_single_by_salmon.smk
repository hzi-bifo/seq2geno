rule gene_counts_single_by_salmon:
    input:
        FQ=lambda wildcards: RNA_READS[wildcards.strain][0],
        REF_SALMON_INDEX_DIR= "{TMP_D}/salmon_index"
    output:
        SALMON_RAW_OUTPUT= directory('{TMP_D}/{strain}/salmon.1')
    params:
        cores=CORES,
        salmon_bin= SOFTWARE['epr_quantifior']
    shell:
        """
        source activate salmon_env
        {params.salmon_bin} quant -i {input.REF_SALMON_INDEX_DIR} \
--gcBias \
-l A -r {input.FQ} \
--sigDigits 0 \
-p {params.cores} \
-o {output.SALMON_RAW_OUTPUT}
        source deactivate
        """

