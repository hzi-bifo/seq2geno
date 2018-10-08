rule vcf_to_indels_per_fam:
    input:
        FAM_VCF=os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}.vcf')       
    output:
        FAM_INDELS_TXT=os.path.join(TMP_D, 'extracted_proteins_nt',
'{fam}_indels.txt'),
        FAM_INDELS_GFF=os.path.join(TMP_D, 'extracted_proteins_nt',
'{fam}_indels.gff'),
        FAM_INDELS_STATS=os.path.join(TMP_D, 'extracted_proteins_nt',
'{fam}_indel_stats.txt'),
    params:
        cores= CORES,
        working_dir='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3',
        strains_perc_cutoff= 0.5,
        len_cutoff= 8
    script: 'vcf2indel.py'

