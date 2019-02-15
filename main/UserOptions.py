'''
Purpose:
    Accept and parse the user options
'''

import os
import argparse


def is_default_file(f):
    outcome= (True if f is None else False)
    return (outcome)

def option_rules(args):
    '''
    Check option conflict and dependency
    '''

def main():

    arg_formatter = lambda prog: argparse.RawTextHelpFormatter(prog,
            max_help_position=4, width = 80)

    parser = argparse.ArgumentParser(
            formatter_class= arg_formatter,
            description='Seq2Geno: the pipline tool '
                'for genomic features computation\n'
                '(Note: all directories and files are relative '
                'to the working directory)')

    parser.add_argument('-v', action= 'version', 
        version='v.Beta')

    ## project
    project_arg= parser.add_argument_group('project')
    project_arg.add_argument('--wd', dest= 'wd', required= True, 
        help='working directory, which every other paths should be relative to')
    project_arg.add_argument('--cores', dest= 'cores', default= 1,
        help='number of cpus')
#    project_arg.add_argument('--mem', dest= 'memory', default= 10,
#        help='max memory size (Gb)')
#    project_arg.add_argument('--log', dest= 'log', default= 'seq2geno.log', 
#        help='working records')
#    project_arg.add_argument('-c', dest='cmprs', action= 'store_true',
#        help='redundancy removal of binary features')
    project_arg.add_argument('--dryrun', dest= 'dryrun', action= 'store_true',
        help='only show the processes')

    ## samples
    sam_arg= parser.add_argument_group('samples')
    sam_arg.add_argument('--dna-reads', dest='dna_reads', type= str,
        help='list of samples and dna sequencing reads', default= '-')
    sam_arg.add_argument('--rna-reads', dest='rna_reads', type= str,
        help='list of samples and rna sequencing reads', default= '-')
    sam_arg.add_argument('--pheno', dest='phe_table', type= str,
        help='list of sample pheno types', default= '-')

    ## reference
    ref_arg= parser.add_argument_group('reference')
    ref_arg.add_argument('--ref-fa', dest='ref_fa', type= str,
        help='reference genome sequences (fasta)', default= '-')
    ref_arg.add_argument('--ref-gbk', dest='ref_gbk', type= str,
        help='reference genome annotation (genbank)', default= '-')
    ref_arg.add_argument('--ref-gff', dest='ref_gff', type= str,
        help='reference genome annotation (gff3)', default= '-')

#    ## snps
#    snps_arg= parser.add_argument_group('snps')
#    snps_arg.add_argument('--snps_f', dest='snps_f', type= str,
#        help='output snps file', default= '-')
#
#    ## denovo
#    denovo_arg= parser.add_argument_group('denovo')
#    denovo_arg.add_argument('--gpa_f', dest='gpa_f', type= str,
#        help='output gene pres/abs file', default= '-')
#
#    ## phylo
#    phylo_arg= parser.add_argument_group('phylo')
#    phylo_arg.add_argument('--tree_f', dest='tree_f', type= str,
#        help='output tree file', default= '-')
#
#    ## expr
#    expr_arg= parser.add_argument_group('expr')
#    expr_arg.add_argument('--expr_f', dest='expr_f', type= str,
#        help='output expression table', default= '-')
#    expr_arg.add_argument('-de', dest='diffexpr', action= 'store_true',
#        help='differential expression analysis')
#    expr_arg.add_argument('-ar', dest='ancrec', action= 'store_true',
#        help='ancestral reconstruction')

    ######
    #####
    args = parser.parse_args()
    args.cores= int(args.cores)
#    print(args)
#    option_rules(args)

    return(args)
