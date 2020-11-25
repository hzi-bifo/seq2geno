#!/usr/bin/env python3
#' Role: Front desk
#' Purpose: 
#' Parse the yaml file and pass the user's options to the package facade
class arguments:
    '''
    The object of arguments, as argparse is replaced
    '''
    def add_opt(self, **entries):
        self.__dict__.update(entries)
    def print_args(self):
        from pprint import pprint
        pprint(vars(self))
    def check_args(self):
        import sys
        # obligatory arguments
        obligatory_args= ['dna_reads', 'wd']
        missing_obligatory_args= [k 
            for k in obligatory_args if not hasattr(self, k)]
        if len(missing_obligatory_args):
            sys.exit('ERROR: "{}" is obligatory'.format(
                ' '.join(missing_obligatory_args)))
        # default values of optional arguments
        optional_args= {'cores':1, 'adaptor': '-', 'rna_reads': '-', 
                'dryrun': True, 
                'phe_table': '-', 'denovo': False, 'snps': False, 'expr': False,
                'phylo': False}
        for k in optional_args:
            if not hasattr(self, k):
                setattr(self, k, optional_args[k])


def parse_arg_yaml(yml_f):
    '''
    Parse the yaml file where the parameters previously were commandline options
    '''
    import yaml
    available_functions= ['snps', 'expr', 'denovo', 'phylo', 'de', 'ar',
    'dryrun']
    with open(yml_f, 'r') as yml_fh:
        opt_dict= yaml.safe_load(yml_fh)
        opt_dict['features']= {k: (True if opt_dict['features'][k] == 'Y' else
                                   False) for k in opt_dict['features']}

    if 'old_config' in opt_dict['general']:
        opt_dict['general']['old_config']= (True if opt_dict['general']['old_config'] == 'Y'
        else False)
    else:
        opt_dict['general']['old_config']= False

    args= arguments()
    try:
        args.add_opt(**opt_dict['general'])
        args.add_opt(**opt_dict['features'])
    except KeyError as e:
        sys.exit('ERROR: {} not found in the input file'.format(str(e)))
    else:
        args.check_args()
        return(args)

def main():
    '''
    Find the yaml file of arguments
    '''
    import yaml 
    import argparse
    import sys

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
    primary_args= parser.parse_args()

    args= parse_arg_yaml(primary_args.yml_f)
    if primary_args.dsply_args:
        args.print_args()
        sys.exit()
    return(args)
