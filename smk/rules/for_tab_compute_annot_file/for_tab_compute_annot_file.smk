rule for_tab_compute_annot_file:
    input:
        ref_gbk=REF_GBK
    output:
        anno_f=temp('annotations_for_snps.tab')
    params:
        species= DNA_READS.index.values.tolist(),
        ref_name= PROJECT
    script:'create_dict_for_snps.py'


