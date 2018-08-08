rule sort_gene_family :
    # compute family-wise alignments
    input:
        roary_gpa= lambda wildcards: os.path.join(TMP_D,
            SOFTWARE['gene_sorter'],'gene_presence_absence.csv'),
        gene_dna_seqs= lambda wildcards: [os.path.join(TMP_D, strain,
            SOFTWARE['assembler'], SOFTWARE['annotator'], 
            'de_novo.ffn') for strain in STRAINS]
    output:
        # slightly different from the original usage; 
        # for the convinience of snakemake
        #output_dir='extracted_proteins_nt/' 
        #seq_files_list=temp('{TMP_D}/{extracted_proteins_dir}/families.list')
        #aln_list='{TMP_D}/aln.list'
        aln= dynamic(os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}.aln'))
    params:
        CORES= CORES,
        MIN_SAMPLES= 2,
        STRAINS= STRAINS,
        TMP_D= TMP_D+'/extracted_proteins_nt'
    script: 'makeGroupAln.py'

'''
rule aln_to_vcf:
    # detect indels from the alignments
    input:
#        core_genes='core_genes_50.txt',
        seq_aln_files_list='{TMP_D}/aln.list'
    output:
        vcf_list= '{TMP_D}/aln_vcf.list'
    params:
        cores= CORES,
        msa2vcf_script='jvarkit/dist/msa2vcf.jar',
        working_dir='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3'
    shell:
#mkdir indels
#cd indels
#        "cd {wildcards.TMP_D}/indels;"
#variant calling on msa
#cat $core_genes | while read i; do echo "java -jar ~/software/jvarkit/dist/msa2vcf.jar< ~/pseudo_genomics/results/assembly/v2/roary/v5/out_95/extracted_proteins_nt/${i}.aln.fasta> ${i}.vcf"; done
#        "cat {params.working_dir}/{input[core_genes]} | while read i;" 
#        "do echo \"java -jar {params.msa2vcf_script} " 
#        "< {params.working_dir}/tmp/extracted_protein_nt/${{i}}.aln.fasta"
#        "> ${{i}}.vcf\";done|"
#        "parallel -j {params.cores};"
#
#        "cd {params.working_dir};"
#        "ls {wildcards.TMP_D}/indels|"
#        "grep vcf$ > {output[vcf_list]} "
        """
        parallel -j {params.cores} \'java -jar {params.msa2vcf_script} \
        < {{}}.aln > {{}}.vcf\' ::: `cat {input[seq_aln_files_list]}| \
        sed \'s/\.aln$//\'`
        cat {input[seq_aln_files_list]} | sed \'s/\.aln$/.vcf/\' \
        > {output[vcf_list]}
        """
'''
rule aln_to_vcf_per_fam:
    input:
        FAM_ALN=os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}.aln')
    output:
        FAM_VCF=os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}.vcf') 
    params:
        cores= CORES,
        msa2vcf_script='jvarkit/dist/msa2vcf.jar',
        working_dir='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3'
    shell:
        """
        java -jar {params.msa2vcf_script} \
        < {input[FAM_ALN]} > {output[FAM_VCF]} 
        """
'''
rule vcf_to_indels:
    input:
        #vcf_list= '{TMP_D}/indels/vcf.list',
        #core_genes='core_genes_50.txt'
        vcf_list= '{TMP_D}/aln_vcf.list'
    output:
        indel_list= '{TMP_D}/indel.list'
    params:
        cores= CORES,
        vcf2indel_script='indel_detection/vcf2indel.py',
        working_dir='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3'
    shell:
#vcf to indel yes/no vector, stats and gff
#cat $core_genes | while read i ; do echo "python /net/sgi/metagenomics/projects/pseudo_genomics/src/PseudoGenomics/indel_detection/vcf2indel.py $i.vcf $i ${i}_indels.txt ${i}_indels.gff ${i}_indel_stats.txt"; done
        """
        source activate py27
        parallel -j {params.cores} \'python {params.vcf2indel_script} \
        {{}}.vcf {{}} {{}}_indels.txt {{}}_indels.gff \
        {{}}_indel_stats.txt \' \
         ::: `cat {input[vcf_list]} | sed \
        \'s/\.vcf$//\'`
        cat {input[vcf_list]} | sed \'s/\.vcf$/.indels.txt/\' \
        > {output[indel_list]}
        source deactivate
        """
'''
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
        vcf2indel_script='indel_detection/vcf2indel.py',
        working_dir='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3'

    shell:
        """
        source activate py27
        python {params.vcf2indel_script} \
        {input} {wildcards.fam} {output[FAM_INDELS_TXT]} \
        {output[FAM_INDELS_GFF]} {output[FAM_INDELS_STATS]} 
        source deactivate
        """

rule expand_by_family:
    input:
        dynamic(os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}_indels.txt'))        
    output:
        indel_list= os.path.join(TMP_D, 'indel.list')
    shell:
        """
        echo \'{input}\' | sed \'s/ /\\n/\' > {output}
        """

rule create_indel_table:
    input:
        roary_gpa_rtab=os.path.join(TMP_D, SOFTWARE['gene_sorter'],
            'gene_presence_absence.Rtab'),
        indel_list= os.path.join(TMP_D, 'indel.list')
    output:
        annot_f=config['indel_table'],
        indel_annot_f='indel_annot.txt',
        indel_stat_f='indels_stats.txt',
        roary_abricate='roary_PA14_abricate.txt'
    params:
        indel_tab_script='indel_detection/generate_indel_features.py'
    shell:
        """
        source activate py27
        python {params.indel_tab_script} \
        <(cut -f1 {input[roary_gpa_rtab]}| tail -n+2 ) \
        {output[annot_f]} {output[indel_annot_f]} \
        {output[indel_stat_f]} {output[roary_abricate]} 
        source deactivate
        """
