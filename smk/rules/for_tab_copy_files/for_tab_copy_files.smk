rule for_tab_copy_files:
    ## duplicate everyone to create 'STRAIN.flt.vcf' and 'STRAIN.flatcount'
    input: 
        flatcount_file=lambda wildcards: os.path.join(TMP_D, wildcards.strain,
'dna_for_tab.flatcount'),
        snp_vcf_file=lambda wildcards: os.path.join(TMP_D, wildcards.strain, 'samtools', 'tab_dna.flt.vcf')
    output:
        flatcount=temp('{strain}.flatcount'),
        flt_vcf=temp('{strain}.flt.vcf')
    wildcard_constraints:
        strain='^[^\/]+'
    shell:
        """
        ln -fs {input.snp_vcf_file} {output.flt_vcf}
        ln -fs {input.flatcount_file} {output.flatcount}
        """
