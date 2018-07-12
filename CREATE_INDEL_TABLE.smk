rule filter_genes: 
    input:
        roary_gpa=TMP_D+'/roary/gene_presence_absence.csv'
    output:
        core_genes='core_genes_50.txt'
    params:
        filter_script='indel_detection/filter.awk' ,
        #awk_threshold=50
        awk_threshold=2
    shell:
    #restrict to genes which are at least in 50 isolate genomes
        "awk -v threshold=\"{params.awk_threshold}\" -f {params.filter_script}  "
        "< {input[roary_gpa]} > {output[core_genes]} "
#    core_genes=~/pseudo_genomics/results/assembly/v2/roary/v5/out_95/core_genes_50.txt

rule sort_gene_family :
    input:
        core_genes='core_genes_50.txt',
        roary_clustered_proteins='{TMP_D}/roary/clustered_proteins',
        gene_dna_seqs= expand('{TMP_D}/prokka/{strain}/{strain}.ffn',
TMP_D=TMP_D, strain= STRAINS)

    output:
        # slightly different from the original usage; 
        # for the convinience of snakemake
        #output_dir='extracted_proteins_nt/' 
        seq_files_list=temp('{TMP_D}/{extracted_proteins_dir}/families.list')
    params:
        make_fasta_script='indel_detection/gene_clusters2multi_fasta.py'
    shell:
        "source activate py27;"
        #number of core genes
        #core_genes_wc=$(wc -l $core_genes | cut -f1 -d" ")
        #extract fasta sequences for each gene family
        #"python {params['sort_gene_family']} {output['output_dir']}"
        "python {params.make_fasta_script} "
        "{wildcards.TMP_D}/{wildcards.extracted_proteins_dir} "
        "<(head -n $(wc -l {input[core_genes]} | cut -f1 -d\" \") "
        "{input[roary_clustered_proteins]}) "
        "{input[gene_dna_seqs]};"
        "ls {wildcards.TMP_D}/{wildcards.extracted_proteins_dir}|"
        "grep fasta$ > {output[seq_files_list]} "

rule make_gene_family_alignment:
    input:
        seq_files_list=temp('{TMP_D}/extracted_protein_nt/families.list')
    output:
        seq_aln_files_list=temp('{TMP_D}/extracted_protein_nt/families_aln.list')
    params:
        cores= CORES,
        working_dir='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3'
    shell:
#cd to roary dir 
#multiple sequence alignment
#mkdir extracted_protein_nt
#cd extracted_protein_nt
        'cd {wildcards.TMP_D}/extracted_protein_nt; '
#for i in `ls`; do i=`echo $i | cut -f1 -d "."` ; echo "mafft $i.fasta > $i.aln.fasta"; done | parallel -j 20
        'parallel -j {params.cores} \'mafft {{}}.fasta > {{}}.aln.fasta\':::'
        '\`ls| grep fasta| cut -f1 -d \".\"\`'
#cd {params.working_dir}/tmp
        'cd {params.working_dir};'
        "ls {wildcards.TMP_D}/indels|"
        "grep aln\.fasta$ > {output[seq_aln_files_list]} "

rule aln_to_vcf:
    input:
        core_genes='core_genes_50.txt',
        seq_aln_files_list='{TMP_D}/extracted_protein_nt/families_aln.list'
    output:
        vcf_list= '{TMP_D}/indels/vcf.list'
    params:
        cores= CORES,
        msa2vcf_script='jvarkit/dist/msa2vcf.jar',
        working_dir='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3'
    shell:
#mkdir indels
#cd indels
        "cd {wildcards.TMP_D}/indels;"
#variant calling on msa
#cat $core_genes | while read i; do echo "java -jar ~/software/jvarkit/dist/msa2vcf.jar< ~/pseudo_genomics/results/assembly/v2/roary/v5/out_95/extracted_proteins_nt/${i}.aln.fasta> ${i}.vcf"; done

        "cat {params.working_dir}/{input[core_genes]} | while read i;" 
        "do echo \"java -jar {params.msa2vcf_script} " 
        "< {params.working_dir}/tmp/extracted_protein_nt/${{i}}.aln.fasta"
        "> ${{i}}.vcf\";done|"
        "parallel -j {params.cores};"

        "cd {params.working_dir};"
        "ls {wildcards.TMP_D}/indels|"
        "grep vcf$ > {output[vcf_list]} "

rule vcf_to_indels:
    input:
        vcf_list= '{TMP_D}/indels/vcf.list',
        core_genes='core_genes_50.txt'
    output:
        indel_list= '{TMP_D}/indels/indel.list'
    params:
        vcf2indel_script='indel_detection/vcf2indel.py',
        working_dir='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3'
    shell:
#vcf to indel yes/no vector, stats and gff
#cat $core_genes | while read i ; do echo "python /net/sgi/metagenomics/projects/pseudo_genomics/src/PseudoGenomics/indel_detection/vcf2indel.py $i.vcf $i ${i}_indels.txt ${i}_indels.gff ${i}_indel_stats.txt"; done
        "cd {wildcards.TMP_D}/indels;"
        "cat {params.working_dir}/{input[core_genes]} | while read i;" 
        "do echo \"java -jar {params.vcf2indel_script} " 
        "< {params.working_dir}/tmp/extracted_protein_nt/${{i}}.aln.fasta"
        "> ${{i}}.vcf\";done;"
        "cd {params.working_dir};"
        "ls {wildcards.TMP_D}/indels|"
        "grep indels.txt$ > {output[indel_list]} "

rule create_indel_table:
    input:
        roary_gpa_rtab=TMP_D+'/roary/gene_presence_absence.Rtab',
        indel_list= TMP_D+'/indels/indel.list'
    output:
        annot_f=config['indel_table'],
        indel_annot_f='indel_annot.txt',
        indel_stat_f='indels_stats.txt',
        roary_abricate='roary_PA14_abricate.txt'
    params:
        indel_tab_script='indel_detection/generate_indel_features.py'
    shell:
        "source activate py27;"
        "python {params.indel_tab_script} <"
        "(cut -f1 {input[roary_gpa_rtab]}| tail -n+2 | grep -v hdl ) "
        "{output[annot_f]} {output[indel_annot_f]} "
        "{output[indel_stat_f]} {output[roary_abricate]}"
