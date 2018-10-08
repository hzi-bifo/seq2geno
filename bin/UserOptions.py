'''
Tasks:
    - Parse the arguments
    - Determine the working folder
    - Check the input files
    - Check the dependencies

Output:
    config file for snakemake
'''

import snakemake
import os
import argparse

def is_default_file(f):
    def_value= '-'
    return (True if f == def_value else False)

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
    functions_arg.add_argument('-w', dest='workdir', type= str, 
        default= os.getcwd(),
        help='set the working directory, where the relative paths are counted')
    functions_arg.add_argument('--use-ng', dest= 'ng', action= 'store_true',
        help='use the faster next-generation version')
    functions_arg.add_argument('-c', dest='comp', action= 'store_true',
        help='compress binary features by pattern')
    functions_arg.add_argument('--dry-run', dest= 'dryrun', action= 'store_true',
        help='dry run the processes')
    functions_arg.add_argument('--keep_temp', dest= 'notemp', action= 'store_true',
        help='do not clean all the intermediate files generated in the process')
    functions_arg.add_argument('--cores', dest= 'cores', default= 1,
        help='number of cpus')

    ## differential expression
    dif_xpr_arg= parser.add_argument_group('differential expression analysis')
    dif_xpr_arg.add_argument('-dx', dest='dif_xpr', action= 'store_true',
        help=
        '''
        detect differentially expressed genes using DESeq2 (Dependency:
        --expr_table, --phe_table)
        ''')
    dif_xpr_arg.add_argument('--dif_xpr_alpha', dest= 'dif_xpr_alpha', 
            default= 0.05, help='the alpha cutoff')
    dif_xpr_arg.add_argument('--dif_xpr_lfc', dest= 'dif_xpr_lfc',
            default= 0, help='the log fold-change cutoff')

    ## ancestral reconstruction of expression levels
    cont_ancrec_arg= parser.add_argument_group('ancestral reconstruction of '\
            'continuous states')
    cont_ancrec_arg.add_argument('-cac', dest='c_ancrec', action= 'store_true',
        help=
        '''
        reconstruct the expression levels along phylogeny using
        phytools::fastAnc (Dependency: --expr_table, --tree)
        ''')

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
    output_arg.add_argument('--tree', nargs= '?', dest='tree',
            default= 'NA',
            type= str, help='the output tree file')
    output_arg.add_argument('--gpa', nargs= '?', dest='gpa_table',
            default= 'NA',
            type= str, help='the output gene pres/abs table')
    output_arg.add_argument('--s-snp', nargs= '?', dest='syn_snps_table',
            default= 'NA',
            type= str, help='the output syn SNPs table')
    output_arg.add_argument('--ns-snp', nargs= '?', dest='nonsyn_snps_table',
            default= 'NA',
            type= str, help='the output non-syn SNPs table')
    output_arg.add_argument('--expr', nargs= '?', dest='expr_table',
            default= 'NA',
            type= str, help='the output expression table')
    output_arg.add_argument('--ind', nargs= '?', dest='indel_table',
            default= 'NA',
            type= str, help='the output indel table')

    ## software
    sw_arg= parser.add_argument_group('software')
    sw_arg.add_argument('--stampy', dest='stampy_exe', type= str,
        help='user-defined version of stampy', 
        default= None)
    sw_arg.add_argument('--raxml', dest='raxml_exe', type= str,
        help='user-defined version of raxml',
        default= None)

    args = parser.parse_args()
    option_rules(args)

    return(args)

'''
## output the config file (yaml format)
config_f= 'test.yaml'
#print(yaml.dump(vars(args), default_flow_style= False))
'''
'''
## create config file
config_txt=''

samples: {}
reference_sequence: {}
reference_annotation: {} 
tree: {}
gpa_table: {}
syn_snps_table: {}
nonsyn_snps_table: {} 
expr_table: {}
indel_table: {}'.format()
'''
'''
w_dir=os.path.abspath('.')
print('working path= {}'.format(w_dir))
seq2geno_dir=os.path.dirname(os.path.realpath(__file__))
print('bin path= {}'.format(seq2geno_dir))
print('config file= {}'.format(args.conf_f))
snakemake.snakemake(snakefile=os.path.join(seq2geno_dir, 'lib/main.smk'),
    configfile= (None if args.conf_f != '' else args.conf_f),
    dryrun=True, workdir= w_dir)


species: 'Paeruginosa'
working_dir: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6'

###### external executives
stampy_exe: '/home/thkuo/bin/stampy-1.0.23/stampy.py'
#raxml_exe: '/home/thkuo/miniconda3/envs/phypal/bin/raxmlHPC-PTHREADS'
#raxml_exe: '/home/thkuo/miniconda3/envs/phypal/bin/raxmlHPC-PTHREADS-AVX2'
raxml_exe: 'raxmlHPC-PTHREADS-AVX2'
#art2genecount_exe: '/net/metagenomics/data/from_moni/old.tzuhao/Paeru/bin/from_SusanneLab/RNA-seq/art2genecount.pl'
#sam2art_exe: '/net/metagenomics/data/from_moni/old.tzuhao/Paeru/bin/from_SusanneLab/RNA-seq/sam2art.pl'

software:
  mapper: 'bwa'
  assembler: 'spades'

cores: 16

##### files
#### sample wise
samples: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/data/samples.tsv'
#### only one
reference_sequence: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/data/reference/RefCln_UCBPP-PA14.fa'
reference_annotation: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/data/reference/RefCln_UCBPP-PA14.gbk'
##roary_gpa:
tree: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/results/Paeru.nwk'
gpa_table: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/results/gpa.tab'
syn_snps_table: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/results/syn_snps.tab'
nonsyn_snps_table: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/results/nonsyn_snps.tab'
expr_table: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/results/expr.tab'
indel_table: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/results/annot.tab'
'''
