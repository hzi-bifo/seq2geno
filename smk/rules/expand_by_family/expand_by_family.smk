rule expand_by_family:
    input:
        FAM_INDELS_FILES=dynamic(os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}_indels.txt'))        
    output:
        indel_list= os.path.join(TMP_D, 'indel.list')
    run:
        out_fh=open(output.indel_list, 'w')
        out_fh.write("\n".join(input.FAM_INDELS_FILES))
        out_fh.close()

