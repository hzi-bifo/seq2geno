'''
Tasks:
    - Parse the arguments
    - Determine the working folder
    - Check the input files
    - Check the dependencies

Output:
    config file for snakemake
'''

import os
import argparse

def is_default_file(f):
    def_value= '-'
    return (True if f is None else False)

def option_rules(args):
    outcome= True
    msg= []
    # needing at least one data
    if all([is_default_file(f) for f in [args.dna_reads, args.rna_reads]]):
        outcome= False
        msg.append('no input data given')

    # differential expression analysis depends on the expression table
    if any(is_default_file(f) for f in [args.dif_xpr, args.expr_table,
        args.phe_table]):
        outcome= False
        msg.append('no expression levels table or phenotype table specified for the differential expression analysis')

    if not outcome:
        exit('\n'.join(msg))

def main():
    parser = argparse.ArgumentParser(
            description='Seq2Geno: the pipline tool '
                'for genomic features computation')

    parser.add_argument('-v', action= 'version', 
        version='v.Beta')

    ## functions
    functions_arg= parser.add_argument_group('functions')
    functions_arg.add_argument('--wd', dest='workdir', type= str, 
        default= os.getcwd(),
        help='set the working directory, where the relative paths are counted')
    functions_arg.add_argument('--ng', dest= 'ng', action= 'store_true',
        help='use the faster next-generation version')
    functions_arg.add_argument('--cmpr', dest='cmpr', action= 'store_true',
        help='compress binary features by pattern')
    functions_arg.add_argument('--dryrun', dest= 'dryrun', action= 'store_true',
        help='dry run the processes')
    functions_arg.add_argument('--keeptemp', dest= 'notemp', action= 'store_true',
        help='do not clean the intermediate files generated in the process')
    functions_arg.add_argument('--cores', dest= 'cores', default= 1,
        help='number of cpus')


    ## software
    sw_arg= parser.add_argument_group('software')
    sw_arg.add_argument('--stampy', dest='stampy_exe', type= str,
        help='user-defined version of stampy', 
        default= None)
    sw_arg.add_argument('--raxml', dest='raxml_exe', type= str,
        help='user-defined version of raxml',
        default= None)


    ## reference
    ref_arg= parser.add_argument_group('reference')
    ref_arg.add_argument('--ref-fa', dest='ref_fa', type= str,
        help='reference genome sequences (fasta)', default= '', required=True)
    ref_arg.add_argument('--ref-gbk', dest='ref_gbk', type= str,
        help='reference genome annotation (genbank)', default= '', required=True)

    ## samples
    sam_arg= parser.add_argument_group('samples')
    sam_arg.add_argument('--dna-reads', dest='dna_reads', type= str,
        help='list of samples and dna sequencing reads', default= '-')
    sam_arg.add_argument('--rna-reads', dest='rna_reads', type= str,
        help='list of samples and rna sequencing reads', default= '-')
    sam_arg.add_argument('--pheno', dest='phe_table', type= str,
        help='list of sample pheno types', default= '-')

    ## outputs
    output_arg= parser.add_argument_group('outputs')
    output_arg.add_argument('--outdir', dest='output_dir',
            default= None,
            type= str, help='the output directory')
    output_arg.add_argument('--tree', dest='tree',
            default= None,
            type= str, 
            help=
            '''
            specifiy a tree filename and compute it if not available
            ''')
    output_arg.add_argument('--gpa', dest='gpa_table',
            default= None,
            type= str, 
            help=
            '''
            specifiy a gene pres/abs table filename and compute it if not
            available
            ''')
    output_arg.add_argument('--s-snp', dest='syn_snps_table',
            default= None,
            type= str, 
            help=
            '''
            specifiy a syn SNPs table filename and compute it if not available
            ''')
    output_arg.add_argument('--ns-snp', dest='nonsyn_snps_table',
            default= None,
            type= str, 
            help=
            '''
            specifiy the non-syn SNPs filename and compute it if not available
            ''')
    output_arg.add_argument('--expr', dest='expr_table',
            default= None,
            type= str, 
            help=
            '''
            specifiy an expression table filename and compute
            it if not available
            ''')
    output_arg.add_argument('--ind', dest='indel_table',
            default= None,
            type= str, 
            help=
            '''
            specifiy an indel table filename and compute it if not available
            ''')

    ## differential expression
    dif_xpr_arg= parser.add_argument_group('differential expression analysis')
    dif_xpr_arg.add_argument('--dx', dest='dif', action= 'store_true',
        help=
        '''
        detect differentially expressed genes using DESeq2 (Dependency:
        --expr_table, --phe_table)
        ''')
    dif_xpr_arg.add_argument('--dif_alpha', dest= 'dif_alpha', 
            default= 0.05, help='the alpha cutoff')
    dif_xpr_arg.add_argument('--dif_lfc', dest= 'dif_lfc',
            default= 0, help='the log fold-change cutoff')
    dif_xpr_arg.add_argument('--dif_out', dest= 'dif_out',
            default= 'difxpr', help='the output folder for analysis results')

    ## ancestral reconstruction of expression levels
    cont_ancrec_arg= parser.add_argument_group('ancestral reconstruction of '\
            'continuous states')
    cont_ancrec_arg.add_argument('--c_ac', dest='c_ac', action= 'store_true',
        help=
        '''
        reconstruct the expression levels along phylogeny using
        phytools::fastAnc (Dependency: --expr_table, --tree)
        ''')
    cont_ancrec_arg.add_argument('--c_ac_out', dest= 'c_ac_out',
            default= 'cont_ancrec', 
            help='the output folder for expression level reconstruction')
  
    ######
    #####
    args = parser.parse_args()
#    option_rules(args)

    return(args)
