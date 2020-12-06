
import unittest
import UserOptions
import seq2geno
from pprint import pprint

def true_parse_arg_yaml(yml_f):
    #' The old method
    '''
    Parse the yaml file where the parameters previously were commandline options
    '''
    import yaml
    available_functions= ['snps', 'expr', 'denovo', 'phylo', 'de', 'ar',
    'dryrun']
    #' read the arguments
    #' flatten the structure
    opt_dict= {}
    with open(yml_f, 'r') as yml_fh:
        opt_dict= yaml.safe_load(yml_fh)
        opt_dict['features']= {k: (True if opt_dict['features'][k] == 'Y' else
                                   False) for k in opt_dict['features']}
    #' reuse the old config files
    if 'old_config' in opt_dict['general']:
        opt_dict['general']['old_config']= (True if opt_dict['general']['old_config'] == 'Y'
        else False)
    else:
        opt_dict['general']['old_config']= False

    args= UserOptions.arguments()
    try:
        args.add_opt(**opt_dict['general'])
        args.add_opt(**opt_dict['features'])
    except KeyError as e:
        sys.exit('ERROR: {} not found in the input file'.format(str(e)))
    else:
        args.check_args()
        return(args)


def true_filter_procs(args):
    import sys
    from tqdm import tqdm
    import os
    from SGProcesses import SGProcess

    config_files= {}
    try:
        ## accept config files
        if args.old_config:
            print('Skip creating config files. '
                  'Old config files will be used if found')
        else:
            import create_config
            config_files= create_config.main(args)
    except:
        print('ERROR: fail to initiate the project')
        e=sys.exc_info()[0]
        print(e)
        sys.exit()

    all_processes= []
    ##>>>
    ## initiate processes

    ## expr
    if args.rna_reads != '-' and args.expr:
        all_processes.append(SGProcess(args.wd,
                      'expr', config_f= config_files['expr'], 
                      dryrun= args.dryrun, 
                      max_cores=
                      int(args.cores)))
    else:
        print('Skip counting expression levels')

    ## snps
    if args.dna_reads != '-' and args.snps:
        all_processes.append(SGProcess(args.wd,
                      'snps', config_f= config_files['snps'], 
                      dryrun= args.dryrun, 
                      max_cores=int(args.cores)))
    else:
        print('Skip calling single nucleotide variants')

    ## denovo
    if args.dna_reads != '-' and args.denovo:
        all_processes.append(SGProcess(args.wd,
                      'denovo', config_f= config_files['denovo'], 
                      dryrun= args.dryrun, 
                      max_cores=int(args.cores)))
    else:
        print('Skip creating de novo assemblies')

    ## phylo
    if args.dna_reads != '-' and args.phylo:
        all_processes.append(SGProcess(args.wd,
                      'phylo', config_f= config_files['phylo'], 
                      dryrun= args.dryrun, 
                      max_cores=int(args.cores)))
    else:
        print('Skip inferring phylogeny')

    ## ancestral reconstruction
    if not args.dryrun and args.dna_reads != '-' and args.rna_reads != '-' and args.ar:
        all_processes.append(SGProcess(args.wd,
                       'ar', config_f= config_files['ar'], 
                       dryrun= args.dryrun, 
                       max_cores=int(args.cores)))
    else:
        print('Skip ancestral reconstruction')

    ## differential expression
    if  not args.dryrun and args.phe_table != '-' and args.rna_reads != '-' and args.de:
        all_processes.append(SGProcess(args.wd,
                      'de', config_f= config_files['de'], 
                       dryrun= args.dryrun, 
                       max_cores=int(args.cores)))
    else:
        print('Skip differential expression analysis')

    return(all_processes)


class TestSelection(unittest.TestCase):
    def test_example_args(self):
        yml_f='seq2geno_inputs.yml'
        #' target to test
        args= UserOptions.parse_arg_yaml(yml_f)
        opted_procs= seq2geno.filter_procs(args)['selected']
        #' true values
        true_args= true_parse_arg_yaml(yml_f)
        true_opted_procs= true_filter_procs(true_args)
        #' comparison
        print('Should be included:')
        print(','.join([p.proc for p in opted_procs]))
        self.assertEqual(set([p.proc for p in opted_procs]),
                set([p.proc for p in true_opted_procs]))

if __name__ == '__main__':
    unittest.main()
