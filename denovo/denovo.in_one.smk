import os
import pandas as pd

list_f=config['list_f']
dna_reads= {}
with open(list_f, 'r') as list_fh:
    for l in list_fh:
        d=l.strip().split('\t')
        dna_reads[d[0]]= d[1].split(',')
out_prokka_dir=config['out_prokka_dir']
out_roary_dir=config['out_roary_dir']
out_spades_dir= config['out_spades_dir']
extracted_proteins_dir=config['extracted_proteins_dir']
out_gpa_f=config['out_gpa_f']
out_indel_f=config['out_indel_f']
REF_GFF=config['REF_GFF']
ref_gbk= config['ref_gbk']
annot_tab=config['annot_tab']
awk_script_f=os.path.join(os.environ['TOOL_HOME'], 'lib', 'filter.awk') 
adaptor_f= config['adaptor']
new_reads_dir= config['new_reads_dir']

rule all:
    input:
        indel_bin_mat= out_indel_f,
        gpa_bin_mat=out_gpa_f,
        non_redundant= expand('{in_tab}_{info}', 
            in_tab= [out_indel_f, out_gpa_f], 
            info= ['GROUPS', 'NONRDNT'])

rule remove_redundant_feat:
    input: 
        F='{in_tab}'
    output: 
        GROUPS='{in_tab}_GROUPS',
        NONRDNT='{in_tab}_NONRDNT'
    conda: 'cmpr_env.yaml'
    script: 'featCompress.py'

rule abricate_dict:
    input:
        roary_clustered_proteins= os.path.join(out_roary_dir,
"clustered_proteins"),
        ref_gff=REF_GFF,
        anno_f=annot_tab
    output:
        rename_dict='{roary_dir}/roary_abricate.txt',
        tmp_annotation_map= '{roary_dir}/annotation_mapped.txt',
        tmp_refined_map= '{roary_dir}/refined_mapping.txt',
        tmp_gff='{roary_dir}/tmp.gff'
    conda:'indel_env.yml'
    shell:
        '''
        cat {input.ref_gff} | sed '/^>/,$d' | tail -n+2 | \
grep -v '#' |grep -v '^\s*$' > {output.tmp_gff}
        field_map_wrapper.edit < <(grep -v @ {input.anno_f}) -s <(cat \
{output.tmp_gff}) -f 6 -m 4 -i| cut -f4,9,19 > {output.tmp_annotation_map}
        refine_mapping.py {output.tmp_annotation_map} \
> {output.tmp_refined_map}
        match_clusters.py {output.tmp_refined_map} \
{input.roary_clustered_proteins} > {output.rename_dict}
        '''

rule gpa_bin_mat:
    input:
        rename_dict_f=os.path.join(out_roary_dir, 'roary_abricate.txt'),
        gpa_csv=os.path.join(out_roary_dir, 'gene_presence_absence.csv')
    output:
        gpa_bin_mat=out_gpa_f
    params:
        strains= list(dna_reads.keys())
    run:
        import pandas as pd
        import textwrap

        strains=params['strains']
        gpa_f=input['gpa_csv']
        output_f=output['gpa_bin_mat']

        ## new name
        name_dict= {}
        in_fh= open(input.rename_dict_f, 'r')
        for l in in_fh:
            d=l.strip().split('\t')
            name_dict[d[0]]= d[1]

        ## read roary result and identify single-copy genes
        df=pd.read_csv(gpa_f, sep= ',',
                header= 0, index_col= 0, quotechar= '"', low_memory=False)

        ## filter and convert the states
        bin_df= df.loc[:, strains].applymap(lambda x: '0' if pd.isna(x) else '1')
        bin_df= bin_df.transpose() # strains in rows
        bin_df.rename(columns=name_dict, inplace= True)
        bin_df.to_csv(output_f, sep= '\t', header= True, index= True,
index_label= 'Gene')
 

rule indel_select_core_genes:
    ## needs shadow
    input:
        gpa_csv=os.path.join('{roary_dir}', 'gene_presence_absence.csv'),
        prot_tab=os.path.join('{roary_dir}', 'clustered_proteins')
    output:
        core_gene_list='{roary_dir}/core_genes_50.txt'
    params:
        min_num_strains= 2, 
        #min_num_strains= 50, 
        filter_awk_script=awk_script_f 
    conda:'indel_env.yml'
    threads:20
    shell:
        '''
        awk -v threshold="{params.min_num_strains}" \
-f {params.filter_awk_script} < {input.gpa_csv} \
| sed 's/\W/_/' \
> {output.core_gene_list} 
        '''

rule indel_align_families:
    input:
        ffn_files=expand(os.path.join(
            out_prokka_dir, '{strain}', '{strain}.ffn'),
            strain=list(dna_reads.keys())),
        core_gene_list=os.path.join(out_roary_dir,'core_genes_50.txt'),
        prot_tab=os.path.join(out_roary_dir, 'clustered_proteins')
    output:
        fam_aln_files=dynamic(os.path.join(
            extracted_proteins_dir, '{fam}.aln'))
    params:
        gene_cluster2multi_script='gene_clusters2multi_fasta.py',
        extracted_proteins_dir=extracted_proteins_dir,
        parallel_log= 'mafft.log'
    threads: 20 
    shell:
        '''
        core_genes={input.core_gene_list}
        #number of core genes
        core_genes_wc=$(wc -l $core_genes | cut -f1 -d" ")
        #extract fasta sequences for each gene family
        {params.gene_cluster2multi_script} \
{params.extracted_proteins_dir} \
<(head -n $core_genes_wc {input.prot_tab}) \
{input.ffn_files}
        # align
        cd {params.extracted_proteins_dir}
        parallel --joblog {params.parallel_log} -j {threads} \
'mafft {{}}.fasta > {{}}.aln' ::: `ls| grep 'fasta'| sed 's/\.fasta//'`
        '''
#for i in `ls|grep \.fasta`; do i=`echo $i | cut -f1 -d "."` ; echo "mafft $i.fasta \
#> $i.aln"; done | parallel --joblog {params.parallel_log} -j {threads}

rule indel_identify_indels:
    input:
        fam_aln_files=dynamic(os.path.join(
            extracted_proteins_dir, '{fam}.aln')),
        core_gene_list=os.path.join(out_roary_dir, 'core_genes_50.txt'),
        gpa_rtab=os.path.join(out_roary_dir, 'gene_presence_absence.Rtab'),
        annot=out_gpa_f,
        roary_abricate= os.path.join(out_roary_dir, 'roary_abricate.txt')
    output:
        indel_annot= out_indel_f
    params:
        vcf2indel_script= 'vcf2indel.py',
        indels_dir='indels',
        generate_feature_script='generate_indel_features.py',
        extracted_proteins_dir=extracted_proteins_dir
    threads: 20
    conda: 'indel_env.yml'
    shell:
        '''
        core_genes=$(cat {input.core_gene_list})
        if [ -d {params.indels_dir} ]; then
            rm -r {params.indels_dir};
        fi 
        mkdir {params.indels_dir} 
        #variant calling on msa
        parallel -j {threads} "msa2vcf \
< {params.extracted_proteins_dir}/{{}}.aln > {params.indels_dir}/{{}}.vcf" ::: $core_genes

        #vcf to indel yes/no vector, stats and gff
        parallel -j {threads} "{params.vcf2indel_script} \
{params.indels_dir}/{{}}.vcf \
{params.indels_dir}/{{}} \
{params.indels_dir}/{{}}_indels.txt \
{params.indels_dir}/{{}}_indels.gff \
{params.indels_dir}/{{}}_indel_stats.txt" ::: $core_genes

        cd {params.indels_dir}
        # In 'clustered_proteins', roary never quotes gene names; 
        # in the .Rtab, space-included names are quoted
        {params.generate_feature_script} \
<(cut -f1 ../{input.gpa_rtab} | tail -n+2 | grep -v hdl | sed 's/"//g') \
../{input.annot} \
../{output.indel_annot} \
../{output.indel_annot}.stats \
../{input.roary_abricate}
        '''

rule roary:
    input:
        gff_files= expand(os.path.join(out_prokka_dir, '{strain}', '{strain}.gff'),
strain= list(dna_reads.keys()))
    output:
        gpa_csv=os.path.join('{roary_dir}', 'gene_presence_absence.csv'),
        gpa_rtab=os.path.join('{roary_dir}', 'gene_presence_absence.Rtab'),
        prot_tab=os.path.join('{roary_dir}', 'clustered_proteins')
    conda: 'perl5_22_env.yml'
    params:
        check_add_perl_env_script= 'install_perl_mods.sh',
        roary_bin= 'roary'
    shadow: "shallow"
    shell:
        '''
        set +u
        {params.check_add_perl_env_script}
        ROARY_HOME=$(dirname $(dirname $(which roary)))
        PERL5LIB=$ROARY_HOME/lib:\
$ROARY_HOME/build/bedtools2/lib:\
$PERL5LIB
        echo $PERL5LIB
        rm -r {wildcards.roary_dir}
        {params.roary_bin} -f {wildcards.roary_dir} \
-e -n -v {input.gff_files} -r -p 30 -g 100000 -z
        set -u
        ''' 

rule create_gff:
    input: os.path.join(out_spades_dir,'{strain}', 'contigs.fasta')
    output: 
        os.path.join(out_prokka_dir, '{strain}', '{strain}.gff'), 
        os.path.join(out_prokka_dir, '{strain}', '{strain}.ffn')
    threads: 1
    conda: 'perl_for_prokka.yml'
    shell:
        '''
        echo $PERL5LIB
        prokka --locustag {wildcards.strain} \
--prefix  {wildcards.strain} \
--force  --cpus {threads} --metagenome --compliant \
--outdir prokka/{wildcards.strain} {input}
        '''
rule create_annot:
    input:
        ref_gbk=ref_gbk
    output:
        anno_f=annot_tab
    params:
        ref_name='reference'
    shell:
        '''
        create_anno.py -r {input.ref_gbk} -n {params.ref_name} -o {output.anno_f}
        '''

rule spades_create_assembly:
    input: 
        READS= lambda wildcards: [
            os.path.join(new_reads_dir,'{}.cleaned.{}.fq.gz'.format(
            wildcards.strain, str(n))) for n in [1,2]]
    output: os.path.join(out_spades_dir,'{strain}', 'contigs.fasta')
    threads:1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000
    params:
        spades_outdir= os.path.join(out_spades_dir, '{strain}'),
        SPADES_OPT='--careful',
        SPADES_BIN='spades.py'
    conda: 'spades_3_10_env.yml'
    script:'run_spades.py'

#rule clean_reads:
#    input:
#        f1= lambda wildcards: os.path.join(
#            new_reads_dir, '{}.fastq.1.gz'.format(wildcards.strain)),
#        f2= lambda wildcards: os.path.join(
#            new_reads_dir, '{}.fastq.2.gz'.format(wildcards.strain))
#    output:
#        log_f= os.path.join(new_reads_dir, '{strain}.log'),
#        f1= os.path.join(new_reads_dir, '{strain}.cleaned.1.fq.gz'),
#        f2= os.path.join(new_reads_dir, '{strain}.cleaned.2.fq.gz')
#    params:
#        adaptor_f= adaptor_f,
#        tmp_f1= lambda wildcards: os.path.join(
#            new_reads_dir, '{}.fastq_cleaned.1'.format(wildcards.strain)),
#        tmp_f2= lambda wildcards: os.path.join(
#            new_reads_dir, '{}.fastq_cleaned.2'.format(wildcards.strain))
#    shadow: "shallow"
#    shell:
#        '''
#        if [ -e "{params.adaptor_f}" ]
#        then
#            fastq-mcf -l 50 -q 20 {params.adaptor_f} {input.f1} {input.f2} \
#-o {params.tmp_f1} -o {params.tmp_f2} > {output.log_f}
#            gzip -9 {params.tmp_f1}
#            gzip -9 {params.tmp_f2}
#        else
#            echo 'No trimming' > {output.log_f}
#            echo $(readlink {input.f1}) >> {output.log_f}
#            echo $(readlink {input.f2}) >> {output.log_f}
#            ln {input.f1} {output.f1}
#            ln {input.f2} {output.f2}
#        fi
#        '''
rule redirect_and_preprocess_reads:
    input: 
        infile1=lambda wildcards: dna_reads[wildcards.strain][0],
        infile2=lambda wildcards: dna_reads[wildcards.strain][1]
    output:
        log_f= os.path.join(new_reads_dir, '{strain}.log'),
        f1= os.path.join(new_reads_dir, '{strain}.cleaned.1.fq.gz'),
        f2= os.path.join(new_reads_dir, '{strain}.cleaned.2.fq.gz')
    params:
        adaptor_f= adaptor_f,
        tmp_f1= lambda wildcards: os.path.join(
            new_reads_dir, '{}.cleaned.1.fq'.format(wildcards.strain)),
        tmp_f2= lambda wildcards: os.path.join(
            new_reads_dir, '{}.cleaned.2.fq'.format(wildcards.strain))
    shell:
         '''
        if [ -e "{params.adaptor_f}" ]
        then
            fastq-mcf -l 50 -q 20 {params.adaptor_f} {input.infile1} {input.infile2} \
-o {params.tmp_f1} -o {params.tmp_f2} > {output.log_f}
            gzip -9 {params.tmp_f1}
            gzip -9 {params.tmp_f2}
        else
            echo 'Reads not trimmed'
            echo 'No trimming' > {output.log_f}
            echo $(readlink {input.infile1}) >> {output.log_f}
            echo $(readlink {input.infile2}) >> {output.log_f}
            ln  -s {input.infile1} {output.f1}
            ln  -s {input.infile2} {output.f2}
        fi
        '''

