## import a list
import pandas as pd
from snakemake.utils import validate

list_f= config['list_f']
#list_f= 'data/samples.10.dna_full.tsv'
dna_reads= {}
with open(list_f, 'r') as list_fh:
    for l in list_fh:
        d=l.strip().split('\t')
        dna_reads[d[0]]= d[1].split(',')

strains= list(dna_reads.keys())

ref_fasta=config['ref_fasta']
ref_gbk=config['ref_gbk']
annot_tab=config['annot_tab']
r_annot=config['r_annot']
snps_table=config['snps_table']
snps_aa_table=config['snps_aa_table']
nonsyn_snps_aa_table=config['nonsyn_snps_aa_table']
snps_aa_bin_mat=config['snps_aa_bin_mat']
nonsyn_snps_aa_bin_mat=config['nonsyn_snps_aa_bin_mat']
adaptor_f= config['adaptor']
new_reads_dir= config['new_reads_dir']

rule all:
    input:
        expand('{strain}.tar.gz', strain= strains),
        snps_aa_bin_mat,
        nonsyn_snps_aa_bin_mat

rule compress_data:
    input:
        ## must ensure the table is done
        snps_aa_bin_mat=snps_aa_bin_mat,
        nonsyn_snps_aa_bin_mat=nonsyn_snps_aa_bin_mat,
        to_compress=['{strain}.sam','{strain}.art','{strain}.sin','{strain}.flatcount',
            '{strain}.rpg','{strain}.stats','{strain}.bam','{strain}.raw.bcf','{strain}.flt.vcf']
    output:
        cmpr_out= '{strain}.tar.gz'
    shell:
        '''
        env GZIP=-9 tar -czvf {output.cmpr_out} {input.to_compress}
        '''

rule create_binary_table:
    input:
        snps_aa_table=snps_aa_table,
        nonsyn_snps_aa_table=nonsyn_snps_aa_table
    output:
        snps_aa_bin_mat=snps_aa_bin_mat,
        nonsyn_snps_aa_bin_mat=nonsyn_snps_aa_bin_mat
    run:
        import pandas as pd
        import re
        def tab_transform(tab_f):
            tab=pd.read_csv(tab_f, sep='\t', header= 0)
            target_columns=['gene', 'pos', 'ref', 'alt', 'ref aa','alt aa']
            tab['feat_name']=tab.apply(lambda r:
    '_'.join(r[target_columns].apply(lambda x: str(x)).tolist()), axis=1)
            mat=tab.set_index('feat_name').drop(target_columns, axis=1).transpose()
            mat=mat.applymap(lambda x: '1' if re.search('[0-9]+', str(x)) else '0')
            return(mat)
        all_mat=tab_transform(input.snps_aa_table)
        all_mat.to_csv(output.snps_aa_bin_mat, sep= '\t', index_label= '')
        nonsyn_mat=tab_transform(input.nonsyn_snps_aa_table)
        nonsyn_mat.to_csv(output.nonsyn_snps_aa_bin_mat, sep= '\t', index_label= '')
        

rule create_table:
    input:
        flt_vcf=expand('{strain}.flt.vcf', strain= strains),
        flatcount=expand('{strain}.flatcount', strain= strains),
        dict_file='dict.txt',
        ref_gbk=ref_gbk,
        annofile=annot_tab
    output:
        snps_table=snps_table
    conda: 'snps_tab_mapping.yml'
    params: 
        tool_script='mutation_table.edit.py'
    shell:
        '''
        {params.tool_script} -f {input.dict_file} -a {input.annofile} -o {output.snps_table}
        '''

rule include_aa_into_table:
    input:
        ref_gbk=ref_gbk,
        snps_table=snps_table
    output:
        snps_aa_table=snps_aa_table,
        nonsyn_snps_aa_table=nonsyn_snps_aa_table
    conda: 'snps_tab_mapping.yml'
    shell:
        '''
        Snp2Amino.py -f {input.snps_table} -g {input.ref_gbk} -o {output.snps_aa_table}
        Snp2Amino.py -n non-syn -f {input.snps_table} -g {input.ref_gbk} \
-o {output.nonsyn_snps_aa_table}
        '''

rule isolate_dict:
    input:
        flt_vcf=expand('{strain}.flt.vcf', strain= strains),
        flatcount=expand('{strain}.flatcount', strain= strains)
    output:
        dict_file='dict.txt'
    threads:1
    wildcard_constraints:
        strain='^[^\/]+$'
    params:
        strains= strains
    run:
        import re
        import os
        ## list and check all required files
        try:
            empty_files= [f for f in input if os.path.getsize(f)==0]
            if len(empty_files) > 0:
                raise Exception('{} should not be empty'.format(
','.join(empty_files)))
        except Exception as e:
            sys.exit(str(e))
        
        with open(output[0], 'w') as out_fh:
            out_fh.write('\n'.join(params.strains))
        
rule my_samtools_SNP_pipeline:
    input:
        sam='{strain}.sam',
        reffile=ref_fasta
    output:
        bam=temp('{strain}.bam'),
        raw_bcf='{strain}.raw.bcf',
        flt_vcf='{strain}.flt.vcf'
    threads:1
    conda: 'snps_tab_mapping.yml'
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/perl5/site_perl/5.22.0:\
$CONDA_PREFIX/lib/perl5/5.22.2:\
$CONDA_PREFIX/lib/perl5/5.22.2/x86_64-linux-thread-multi/:\
$PERL5LIB
        echo $PERL5LIB
        my_samtools_SNP_pipeline {wildcards.strain} {input.reffile} 0
        """

rule my_stampy_pipeline_PE:
    input:
        infile1= lambda wildcards: os.path.join(
        new_reads_dir, '{}.fastq_cleaned.1.gz'.format(wildcards.strain)),
        infile2= lambda wildcards: os.path.join(
        new_reads_dir, '{}.fastq_cleaned.2.gz'.format(wildcards.strain)),
        reffile=ref_fasta,
        ref_index_stampy=ref_fasta+'.stidx',
        ref_index_bwa=ref_fasta+'.bwt',
        annofile=annot_tab,
        Rannofile=r_annot
    output:
        sam=temp('{strain}.sam'),
        art='{strain}.art',
        sin='{strain}.sin',
        flatcount='{strain}.flatcount',
        rpg='{strain}.rpg',
        stat='{strain}.stats'
    threads:1
    conda: 'snps_tab_mapping.yml'
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/perl5/5.22.2/x86_64-linux-thread-multi/:$PERL5LIB
        export PERL5LIB=$CONDA_PREFIX/lib/perl5/5.22.2:$PERL5LIB
        export PERL5LIB=$CONDA_PREFIX/lib/perl5/site_perl/5.22.0:$PERL5LIB
        my_stampy_pipeline_PE {wildcards.strain} \
{input.infile1} {input.infile2} {input.reffile} \
{input.annofile} {input.Rannofile} 2> {wildcards.strain}.log
        """

#rule clean_reads:
#    input:
#        f1= lambda wildcards: os.path.join(
#            new_reads_dir, '{}.fastq.1.gz'.format(wildcards.strain)),
#        f2= lambda wildcards: os.path.join(
#            new_reads_dir, '{}.fastq.2.gz'.format(wildcards.strain))
#    output:
#        log_f= os.path.join(new_reads_dir, '{strain}.log'),
#        f1= os.path.join(new_reads_dir, '{strain}.fastq_cleaned.1.gz'),
#        f2= os.path.join(new_reads_dir, '{strain}.fastq_cleaned.2.gz')
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
        f1= os.path.join(new_reads_dir, '{strain}.fastq_cleaned.1.gz'),
        f2= os.path.join(new_reads_dir, '{strain}.fastq_cleaned.2.gz')
    params:
        adaptor_f= adaptor_f,
        tmp_f1= lambda wildcards: os.path.join(
            new_reads_dir, '{}.fastq_cleaned.1'.format(wildcards.strain)),
        tmp_f2= lambda wildcards: os.path.join(
            new_reads_dir, '{}.fastq_cleaned.2'.format(wildcards.strain))
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

rule create_r_annot:
    input:
        ref_gbk=ref_gbk
    output:
        R_anno_f=r_annot
    shell:
        '''
        create_R_anno.py -r {input.ref_gbk} -o {output.R_anno_f}
        '''

rule stampy_index_ref:
    input:
        reffile=ref_fasta
    output:
        ref_fasta+'.bwt',
        ref_fasta+'.stidx',
        ref_fasta+'.sthash'
    conda: 'snps_tab_mapping.yml'
    shell:
        '''
        stampy.py -G {input.reffile} {input.reffile}
        stampy.py -g {input.reffile} -H {input.reffile}
        bwa index -a bwtsw {input.reffile}
        '''
