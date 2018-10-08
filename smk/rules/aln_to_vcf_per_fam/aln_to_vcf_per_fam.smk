rule aln_to_vcf_per_fam:
    input:
        FAM_ALN=os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}.aln')
    output:
        FAM_VCF=os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}.vcf') 
    params:
        cores= CORES,
        msa2vcf_script='lib/jvarkit/dist/msa2vcf.jar',
        working_dir='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3'
    shell:
        """
        java -jar {params.msa2vcf_script} \
        < {input[FAM_ALN]} > {output[FAM_VCF]} 
        """

