rule for_tab_compute_annot_file:
    input:
        ref_gbk=config['reference_annotation']
    output:
        anno_f=temp('annotations_for_snps.tab')
    params:
        species= config['species']
    script:'create_dict_for_snps.py'


