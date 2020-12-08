#!/usr/bin/env python3
#' Role: Front desk
#' Purpose: 
#' Parse the yaml file and pass the user's options to the package facade
import yaml
import os
from pprint import pprint
import sys
import argparse
import sys
class arguments:
    '''
    The object of arguments, as argparse is replaced
    '''
    def add_opt(self, **entries):
        self.__dict__.update(entries)
    def print_args(self):
        pprint(vars(self))
    def check_args(self):
        import ArgsTest

        #' default values of optional arguments
        optional_args= {'cores':1, 'mem_mb':-1, 
                        'adaptor': '-', 'rna_reads': '-', 
                'dryrun': 'Y', 'phe_table': '-', 
                'denovo': 'N', 'snps': 'N', 'expr': 'N', 'phylo': 'N'}
        for k in optional_args:
            if not hasattr(self, k):
                #' blanck
                setattr(self, k, optional_args[k])
            elif len(str(getattr(self, k))) == 0:
                #' empty 
                setattr(self, k, optional_args[k])


        #' obligatory arguments
        obligatory_args= ['dna_reads', 'ref_fa', 'wd']
        #' ensure obligatory arguments included
        for  k in obligatory_args:
            assert hasattr(self, k), 'Obligatory arguments "{}" not properly set'.format(
                    k)
        #' check the reference genome
        ArgsTest.test_reference_seq(self.ref_fa)
        #' check the dna-seq data
        ArgsTest.test_dna_reads_list(self.dna_reads)
        #' check the rna-seq data (if set)
        if hasattr(self, 'rna_reads'):
            ArgsTest.test_rna_reads_list(self.rna_reads)



def parse_arg_yaml(yml_f):
    '''
    Parse the yaml file where the parameters previously were commandline options
    '''
    available_functions= ['snps', 'expr', 'denovo', 'phylo', 'de', 'ar',
    'dryrun']
    #' read the arguments
    #' flatten the structure
    opt_dict= {}
    with open(yml_f, 'r') as yml_fh:
        opt_dict= yaml.safe_load(yml_fh)
    #' reuse the old config files
    if not ('old_config' in opt_dict['general']):
        opt_dict['general']['old_config']= 'N'

    args= arguments()
    try:
        args.add_opt(**opt_dict['general'])
        args.add_opt(**opt_dict['features'])
    except KeyError as e:
        sys.exit('ERROR: {} not found in the input file'.format(str(e)))
    else:
        args.check_args()
        return(args)

def check_primary_args(primary_args):
    #' check the primary argument
    #' the log file must not exist
    log_f= primary_args.log_f
    if not primary_args.log_f is None:
        # The log file should either be skipped, which will have the messages
        # directed to stdout and stderr as usual, or be a non-used filename,
        # where the stdout and stderr will be redirected
        assert (not os.path.isfile(getattr(primary_args, 'log_f'))), 'Log file existing. Please specify another filename'
        print('Log recorded in {}'.format(log_f))
        sys.stdout = open(log_f, 'w')
        sys.stderr = sys.stdout 

    #' the yml file must exist
    assert os.path.isfile(primary_args.yml_f), 'The yaml file not existing'
    print('#CONFIGFILE:{}'.format(primary_args.yml_f))

def main():
    '''
    Find the yaml file of arguments
    '''

    arg_formatter = lambda prog: argparse.RawTextHelpFormatter(prog,
            max_help_position=4, width = 80)

    parser = argparse.ArgumentParser(
            formatter_class= arg_formatter,
            description='''
        Seq2Geno: the automatic tool for computing genomic features from
        the sequencing data\n''')

    parser.add_argument('-v', action= 'version', 
        version='v.Beta')
    parser.add_argument('-d', dest= 'dsply_args', action= 'store_true',
        help= 'show the arguments described in the config file (yaml) and exit')
    parser.add_argument('-f', dest= 'yml_f', required= True, 
        help= 'the yaml file where the arguments are listed')
    parser.add_argument('-l', dest= 'log_f', required= False, 
        help= 'a non-existing filename for log')
    primary_args= parser.parse_args()

    #' check those primary arguments
    check_primary_args(primary_args)

    args= parse_arg_yaml(primary_args.yml_f)
    #' display the primary arguments only 
    if primary_args.dsply_args:
        args.print_args()
        sys.exit(0)

    return(args)
