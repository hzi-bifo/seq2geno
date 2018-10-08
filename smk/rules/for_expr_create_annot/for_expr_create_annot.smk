rule for_expr_create_annot:
    input:
        ref_gbk=config['ref_gbk']
    output:
        R_anno_f=temp('{TMP_D}/R_annotations.tab')
    script:'create_R_anno.py'

