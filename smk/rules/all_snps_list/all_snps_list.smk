rule all_snps_list:
    input:
        ## the exact files should also be, so they cannot be deleted before 
        flatcount_file=expand(os.path.join(TMP_D, '{strain}', 'dna_for_tab.flatcount'),
strain=DNA_READS.index.values.tolist()),
        snp_vcf_file=expand(os.path.join(TMP_D, '{strain}', 'samtools',
'tab_dna.flt.vcf'), strain=DNA_READS.index.values.tolist()),
        flt_files=expand("{strain}.flatcount", 
            strain=DNA_READS.index.values.tolist()),
        flt_vcf_files=expand("{strain}.flt.vcf", 
            strain=DNA_READS.index.values.tolist()),
        dict_f= os.path.join(TMP_D, 'dict.txt'),
        anno_f='annotations_for_snps.tab'
    output: 
        snps_list=temp(os.path.join(TMP_D, 'all_SNPs.tab'))

    conda: ENV_FILES_POOL.find_yaml('old_mapping')
    threads: 1
    params:
        script_f='mutation_table.py'
    shell:
        """
        {params.script_f} -f {input.dict_f} \
-a {input.anno_f} \
-o {output.snps_list}
        """

