include: os.path.join(RULE_LIB_DIR, 'for_tab_filter_vcf', 'for_tab_filter_vcf.smk')
include: os.path.join(RULE_LIB_DIR, 'for_tab_create_vcf', 'for_tab_create_vcf.smk')
include: os.path.join(RULE_LIB_DIR, 'for_tab_sam2bam', 'for_tab_sam2bam.smk')
include: os.path.join(RULE_LIB_DIR, 'for_snps_convert_to_art', 'for_snps_convert_to_art.smk')
include: os.path.join(RULE_LIB_DIR, 'for_snps_convert_to_sin', 'for_snps_convert_to_sin.smk')
include: os.path.join(RULE_LIB_DIR, 'for_snps_convert_to_flatcount', 'for_snps_convert_to_flatcount.smk')
include: os.path.join(RULE_LIB_DIR, 'for_tab_redirect_stampy_result', 'for_tab_redirect_stampy_result.smk')
include: os.path.join(RULE_LIB_DIR, 'for_tab_stampy_single_mapping', 'for_tab_stampy_single_mapping.smk')
include: os.path.join(RULE_LIB_DIR, 'for_tab_stampy_paired_mapping', 'for_tab_stampy_paired_mapping.smk')