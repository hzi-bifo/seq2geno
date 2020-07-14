#!/usr/bin/env python3

'''
Purpose:
    Interact with the user, accept the request, and convert the information
    to the input for the subsequent snakemake workflow
'''

import os
import subprocess
import sys


class phylo_args:
    def __init__(self, list_f, REF_FA, REF_GFF,adaptor,
        redirected_reads_dir, config_f='config.yml'):
        self.list_f= list_f
        self.REF_FA= REF_FA
        self.REF_GFF= REF_GFF
        self.adaptor= adaptor
        self.new_reads_dir= redirected_reads_dir
        subfolder= 'phylo'
        self.config_f= os.path.join(
            subfolder, config_f)
        self.seq_list= 'seq_list'
        self.results_dir= 'tmp'
        self.families_seq_dir= 'families'
        self.aln_f= 'OneLarge.gapReplaced.var2.gt90.aln'
        self.tree_f= 'OneLarge.gapReplaced.var2.gt90.nwk'

class expr_args:
    def __init__(self, list_f, ref_fasta, ref_gbk,  
    config_f='config.yml'):
        self.list_f= list_f
        self.ref_fasta= ref_fasta
        self.ref_gbk= ref_gbk 
        #self.out_table=os.path.basename(out_table)
        #self.out_log_table= self.out_table+'_log_transformed'
        #self.config_f= os.path.join(os.path.dirname(out_table), 
        #        config_f)
        self.out_table='expr.mat'
        self.out_log_table= 'expr.mat_log_transformed'
        subfolder= 'expr'
        self.config_f= os.path.join(
            subfolder, config_f)
        #self.annot_tab='data/Pseudomonas_aeruginosa_PA14_annotation_with_ncRNAs_07_2011_12genes.tab'
        #self.r_annot= 'data/Pseudomonas_aeruginosa_PA14_12genes_R_annotation'
        self.annot_tab='annotation.tab'
        self.r_annot='R_annotation.tab'

class snps_args:
    def __init__(self, list_f, ref_fasta, ref_gbk, adaptor,
            redirected_reads_dir, config_f= 'config.yml'):
        self.list_f=list_f
        self.adaptor= adaptor
        self.new_reads_dir= redirected_reads_dir
        self.ref_fasta= ref_fasta
        self.ref_gbk= ref_gbk 
        #self.snps_aa_bin_mat=snps_aa_bin_mat
        #self.nonsyn_snps_aa_bin_mat= 'nonsyn_'+ self.snps_aa_bin_mat
        #self.config_f= os.path.join(
        #    os.path.dirname(snps_aa_bin_mat), config_f)
        subfolder= 'snps'
        self.config_f= os.path.join(
            subfolder, config_f)
        #self.annot_tab= 'data/Pseudomonas_aeruginosa_PA14_annotation_with_ncRNAs_07_2011_12genes.tab'
        #self.r_annot= 'data/Pseudomonas_aeruginosa_PA14_12genes_R_annotation'
        #self.snps_table= 'all_SNPs.tab'
        #self.snps_aa_table= 'all_SNPs_final.tab'
        #self.nonsyn_snps_aa_table= 'nonsyn_SNPs_final.tab'
        self.annot_tab='annotation.tab'
        self.r_annot='R_annotation.tab'
        self.snps_table='all_SNPs.tab'
        self.snps_aa_table='all_SNPs_final.tab'
        self.nonsyn_snps_aa_table='nonsyn_SNPs_final.tab'
        self.snps_aa_bin_mat='all_SNPs_final.bin.mat'
        self.nonsyn_snps_aa_bin_mat='nonsyn_SNPs_final.bin.mat'

class denovo_args:
    def __init__(self, list_f, REF_GFF, ref_gbk, adaptor,
        redirected_reads_dir, config_f='config.yml'):
        self.list_f= list_f
        self.REF_GFF= REF_GFF
        self.ref_gbk= ref_gbk 
        self.new_reads_dir= redirected_reads_dir
        self.adaptor= adaptor
        self.out_spades_dir= 'spades'
        self.out_prokka_dir= 'prokka'
        self.out_roary_dir= 'roary'
        self.extracted_proteins_dir= 'extracted_proteins_nt/'
        #self.out_gpa_f= os.path.basename(out_gpa_f)
        #self.out_indel_f= 'indel_'+self.out_gpa_f
        self.out_gpa_f= 'gpa.mat'
        self.out_indel_f= 'indel.mat'
        #self.anno_f= 'data/Pseudomonas_aeruginosa_PA14_annotation_with_ncRNAs_07_2011_12genes.tab'
        self.annot_tab='annotation.tab'
        subfolder= 'denovo'
        self.config_f= os.path.join(
            subfolder, config_f)

class cmprs_args:
    def __init__(self, bin_tables,
        config_f='config.yml'):
        self.config_f= config_f
        self.in_tab=','.join(bin_tables)

class ar_args:
    def __init__(self, phylo_config_f, expr_config_f,
        config_f='config.yml'):
        self.config_f= config_f
        self.phylo_config= phylo_config_f
        self.expr_config= expr_config_f
        self.C_ANCREC_OUT= os.path.join(os.path.dirname(expr_config_f), 
            'ancrec')

class diffexpr_args:
    def __init__(self, pheno, expr_config_f,
        config_f='config.yml'):
        self.config_f= config_f
        self.pheno_tab= pheno
        self.expr_config= expr_config_f
        self.out_dir= os.path.join(os.path.dirname(expr_config_f), 'dif')
        self.alpha= 0.05
        self.lfc= 1

def create_yaml_f(args, wd, config_f):
    import yaml
    import os
    try:
        config_f= os.path.abspath(os.path.join(wd, config_f))
        print(config_f)
        target_dir= os.path.dirname(config_f)
        if not os.path.exists(target_dir):
            print('creating {}...'.format(target_dir))
            os.makedirs(target_dir)
        config_fh= open(config_f, 'w')
        yaml.dump(vars(args), config_fh, default_flow_style= False)
        config_fh.close()
    except Exception as e:
        from datetime import datetime
        print('ERROR ({})'.format(proc))
        print('{}\t{}\n'.format(
            datetime.now().isoformat(' ',timespec= 'minutes'),
            str(e)))
        sys.exit()



def main(args):
    print('Creating config files')
#    seq2geno_home= os.path.abspath(
#            os.path.join(os.path.dirname(os.path.realpath(__file__)),
#                os.pardir))
#
#    # parse user's arguments
#    import UserOptions
#    args= UserOptions.main()
    # max core number
    cores= args.cores

    # snps
    s_args= snps_args(list_f= os.path.abspath(args.dna_reads), 
            adaptor= args.adaptor,
            redirected_reads_dir= os.path.join(
                os.path.abspath(args.wd), 'reads', 'dna'),
            ref_fasta= args.ref_fa, 
            ref_gbk= args.ref_gbk)
    create_yaml_f(s_args, args.wd, s_args.config_f)

    # denovo
    d_args= denovo_args(list_f= os.path.abspath(args.dna_reads), 
            REF_GFF= args.ref_gff, 
            adaptor= args.adaptor,
            redirected_reads_dir= os.path.join(
                os.path.abspath(args.wd), 'reads', 'dna'),
            ref_gbk= args.ref_gbk)
    create_yaml_f(d_args, args.wd, d_args.config_f)

    # phylo
    p_args= phylo_args(list_f= os.path.abspath(args.dna_reads), 
            adaptor= args.adaptor,
            redirected_reads_dir= os.path.join(
                os.path.abspath(args.wd), 'reads', 'dna'),
            REF_FA= args.ref_fa, 
            REF_GFF= args.ref_gff)
    create_yaml_f(p_args, args.wd, p_args.config_f)

    # expr
    e_args= expr_args(list_f= os.path.abspath(args.rna_reads), 
            ref_fasta= args.ref_fa, 
            ref_gbk= args.ref_gbk)
    create_yaml_f(e_args, args.wd, e_args.config_f)

    # compress
    c_args= cmprs_args(bin_tables=[
        os.path.join(os.path.dirname(s_args.config_f),
            s_args.snps_aa_bin_mat),
        os.path.join(os.path.dirname(s_args.config_f),
            s_args.nonsyn_snps_aa_bin_mat), 
        os.path.join(os.path.dirname(d_args.config_f), 
            d_args.out_gpa_f),
        os.path.join(os.path.dirname(d_args.config_f), 
            d_args.out_indel_f)])
    create_yaml_f(c_args, args.wd, c_args.config_f)

    # ancestral reconstruction
    a_args= ar_args(phylo_config_f= p_args.config_f,
            expr_config_f= e_args.config_f)
    create_yaml_f(a_args, args.wd, a_args.config_f)
    
    # differential expression analysis
    de_args= diffexpr_args(pheno= os.path.abspath(args.phe_table),
            expr_config_f= e_args.config_f)
    create_yaml_f(de_args, args.wd, de_args.config_f)

    # Determine which version to use (ori or ng)
    # Create the environment, followed by the analysis snakemake workflows
    #subprocess.call(main_cmd)
    return({'snps':  os.path.abspath(os.path.join(args.wd, s_args.config_f)),
            'expr': os.path.abspath(os.path.join(args.wd, e_args.config_f)),
            'denovo':os.path.abspath(os.path.join(args.wd, d_args.config_f)),
            'phylo':os.path.abspath(os.path.join(args.wd, p_args.config_f)),
            'ar':os.path.abspath(os.path.join(args.wd, a_args.config_f)),
            'de':os.path.abspath(os.path.join(args.wd, de_args.config_f))})
