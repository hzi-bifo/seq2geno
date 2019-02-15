rule compute_gpa_table:
    input:
        rename_dict=TMP_D+'/prokkaPA14.txt',
        roary_gpa=TMP_D+"/roary/gene_presence_absence.csv"
    output:
        gpa_table=GPA_OUT
    params: 
        strains= DNA_READS.index.values.tolist()
    script: 'roary_gpa2bin.py'

rule compute_rename_dict:
    input:
        roary_clustered_proteins= TMP_D+"/roary/clustered_proteins",
        ref_gff=REF_GFF,
        anno_f='annotations_for_snps.tab' 
    output:
        rename_dict=temp(TMP_D+'/prokkaPA14.txt'),
        tmp_annotation_map= temp(TMP_D+'/annotation_mapped.txt'),
        tmp_refined_map= temp(TMP_D+'/refined_mapping.txt'),
        tmp_gff= temp(TMP_D+'/tmp.gff')
    conda:ENV_FILES_POOL.find_yaml('gpa_rename_env')
    shell:
        '''
        cat {input.ref_gff} | sed '/^>/,$d' | tail -n+2 | \
grep -v '#' > {output.tmp_gff}
        ./field_map_wrapper.edit < <(grep -v @ {input.annot_f}) -s <(cat \
{output.tmp_gff}) -f 6 -m 4 -i| cut -f4,9,19 > {output.tmp_annotation_map}
        python refine_mapping.py {poutput.tmp_annotation_map} \
> {output.tmp_refined_map}
        python match_clusters.py {output.tmp_refined_map} \
{input.roary_clustered_proteins} > {output.rename_dict}
        '''
