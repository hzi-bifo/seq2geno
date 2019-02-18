#!/usr/bin/env python3

def collect_results(project_dir):
    import shutil
    import yaml
    results_newdir= os.path.join(project_dir, 'RESULTS')
    if os.path.isdir(results_newdir):
        shutil.rmtree(results_newdir)
    os.mkdir(results_newdir)
    denovo_config= os.path.join(project_dir, 'denovo', 'denovo_config.yml') 
    snps_config= os.path.join(project_dir, 'snps', 'snps_config.yml') 
    expr_config= os.path.join(project_dir, 'expr', 'expr_config.yml') 
    phylo_config= os.path.join(project_dir, 'phylo', 'phylo_config.yml') 
    diffexpr_config= os.path.join(project_dir, 'de_config.yml')

    ## de novo assemblies
    d_config= yaml.safe_load(open(denovo_config, 'r'))
    assem_newdir= os.path.join(results_newdir, 'assemblies')
    os.mkdir(assem_newdir)
    assem_outdir= os.path.join(os.path.dirname(denovo_config),
            d_config['out_spades_dir'])
    strains=os.listdir(assem_outdir)
    for s in strains:
        old= os.path.join(assem_outdir, s, 'contigs.fasta')
        new= os.path.join(assem_newdir, '{}.fasta'.format(s))
        os.symlink(os.path.abspath(old), new)
    ## gpa and indel
    bin_tab_newdir= os.path.join(results_newdir, 'bin_tables')
    os.mkdir(bin_tab_newdir)
    nr_gpa_mat=os.path.join(os.path.dirname(denovo_config), 
        d_config['out_gpa_f'])+'_NONRDNT'
    nr_indel_mat= os.path.join(os.path.dirname(denovo_config), 
        d_config['out_indel_f'])+'_NONRDNT' 
    ## snps
    s_config= yaml.safe_load(open(snps_config, 'r'))
    nr_all_snps_mat= os.path.join(os.path.dirname(snps_config), 
        s_config['snps_aa_bin_mat'])+'_NONRDNT'
    nr_ns_snps_mat= os.path.join(os.path.dirname(snps_config), 
        s_config['nonsyn_snps_aa_bin_mat'])+'_NONRDNT'
    for f in [nr_gpa_mat, nr_indel_mat, nr_all_snps_mat, nr_ns_snps_mat]:
        new= os.path.join(bin_tab_newdir, os.path.basename(f))
        os.symlink(os.path.abspath(f), new)
    ## expr
    num_tab_newdir= os.path.join(results_newdir, 'num_tables')
    os.mkdir(num_tab_newdir)
    e_config= yaml.safe_load(open(expr_config, 'r'))
    expr_mat= os.path.join(os.path.dirname(expr_config), 
        e_config['out_log_table'])
    new= os.path.join(num_tab_newdir, 'expr.log.mat')
    os.symlink(os.path.abspath(expr_mat), new)

    ar_outdir= os.path.join(os.path.dirname(expr_config), 'ancrec')
    ar_newdir= os.path.join(results_newdir, 'expr_ancestral_reconstructions')
    os.symlink(os.path.abspath(ar_outdir), ar_newdir)

    de_outdir= os.path.join(os.path.dirname(expr_config), 'dif/')
    de_newdir= os.path.join(results_newdir, 'differential_expression')
    os.symlink(os.path.abspath(de_outdir), de_newdir)
    
    ## phylogeny
    p_config= yaml.safe_load(open(phylo_config, 'r'))
    phy_newdir= os.path.join(results_newdir, 'phylogeny')
    os.mkdir(phy_newdir)
    phylo_nwk= os.path.join(os.path.dirname(phylo_config), 
        p_config['tree_f'])
    new= os.path.join(phy_newdir, 'tree.nwk')
    os.symlink(os.path.abspath(phylo_nwk), new)
    ## phenotype table
    de_config= yaml.safe_load(open(diffexpr_config, 'r'))
    phe_mat= de_config['pheno_tab']
    phe_newdir= os.path.join(results_newdir, 'phenotype')
    os.mkdir(phe_newdir)
    new= os.path.join(phe_newdir, 'phenotypes.mat')
    os.symlink(os.path.abspath(os.path.realpath(phe_mat)), new)


def EditEnv(proc):
    import re
    import pandas as pd
    from datetime import datetime
    script_dir=os.path.dirname(os.path.realpath(__file__))
    toolpaths_f=os.path.join(script_dir, 'ToolPaths.tsv')
    ## read the env variables
    env_df= pd.read_csv(toolpaths_f, sep= '\t', comment= '#', index_col= 0)
    env_series=pd.Series([])
    env_dict= {}
    try:
        env_series= env_df.loc[proc,:]
    except KeyError as ke:
        print('ERROR ({})'.format(proc))
        print('{}\t{}\n'.format(
            datetime.now().isoformat(' ',timespec= 'minutes'),
            'unavailable function'))
        sys.exit()
    else:
        try:
            os.environ['TOOL_HOME']= os.path.join(os.environ['SEQ2GENO_HOME'],
                    str(env_series['TOOL_HOME']))
            all_env_var= dict(os.environ)
            env_dict={'TOOL_HOME': all_env_var['TOOL_HOME']}
            for env in env_series.index.values.tolist():
                if env == 'TOOL_HOME':
                    continue
                val= str(env_series[env])
                included= list(set(re.findall('\$(\w+)', val)))
                for included_env in included:
                    val= re.sub('\$'+included_env, all_env_var[included_env], val)
                env_dict[env]= val
        except:
            print('ERROR ({})'.format(proc))
            print('{}\t{}\n'.format(
                datetime.now().isoformat(' ',timespec= 'minutes'),
                'Unable to set environment variables'))
            sys.exit()
    return(env_dict)

def run_proc(proc, config_f, dryrun= True, max_cores=1):
    env_dict=EditEnv(proc)
    print(proc)
    try:
        import snakemake
        import os 
        os.environ['PATH']=env_dict['PATH']
        success= snakemake.snakemake(
            configfile=config_f,
            snakefile= env_dict['SNAKEFILE'],
            workdir= os.path.dirname(config_f),
            unlock= True)
        if not success:
            raise Exception('Could not unlock target folder'.format(os.path.dirname(config_f)))

        success=snakemake.snakemake(
            snakefile= env_dict['SNAKEFILE'],
            cores= max_cores,
            configfile=config_f,
            workdir= os.path.dirname(config_f),
            use_conda=True,
            conda_prefix= os.path.join(env_dict['TOOL_HOME'], 'env'),
            dryrun= dryrun,
            printshellcmds= True,
            force_incomplete= True,
            notemp=True
            )
        if not success:
            raise Exception('Snakemake workflow fails')
    except Exception as e:
        from datetime import datetime
        print('ERROR ({})'.format(proc))
        print('{}\t{}\n'.format(
            datetime.now().isoformat(' ',timespec= 'minutes'),
            e))
        raise RuntimeError(e)
    except :
        from datetime import datetime
        print('ERROR ({})'.format(proc))
        print('{}\t{}\n'.format(
            datetime.now().isoformat(' ',timespec= 'minutes'),
            sys.exc_info()))
        raise RuntimeError('Unknown problem occured when lauching Snakemake')

#proc_targets= {'denovo': 'denovo.log', 'snps': 'snps.log'}
#config_f= {'denovo': 'denovo/config.yml'}
#targets= []
#required= ['denovo', 'snps']
#max_cores= 15
#dryrun= True
#all_analysis= ['denovo', 'snps', 'phylo', 'expr']

from pprint import pprint
import os
import sys

import UserOptions
args= UserOptions.main()
from pprint import pprint
pprint(args)
try:
    import create_config
    create_config.main(args)
except:
    print('ERROR: fail to initiate the project')
    sys.exit()

'''
import argparse
arg_formatter = lambda prog: argparse.RawTextHelpFormatter(prog,
        max_help_position=4, width = 80)
parser = argparse.ArgumentParser(
        formatter_class= arg_formatter,
        description='start seq2geno')

parser.add_argument('-v', action= 'version', 
    version='v.Beta')
parser.add_argument('-f', dest= 'target_dir', type= str,
    help= 'target folder, created by create_config', 
    required= True)
parser.add_argument('-c', dest= 'max_core', type= int,
    help= 'max number of cores to use', 
    default= 1)
#parser.add_argument('-m', dest= 'max_mem_per_task', type= int,
#    help= 'max size (Gb) of memory for each task', 
#    default= 10)
parser.add_argument('--dryrun', dest= 'dryrun', action= 'store_true',
    help= 'only list the computational processes', default= False)
args= parser.parse_args()
'''
try:
    ## expr
    if args.rna_reads != '-':
        config_f= os.path.join(args.wd, 'expr', 'expr_config.yml')
        run_proc('expr', config_f, dryrun= args.dryrun, max_cores= args.cores)
    else:
        print('Expression level will not be counted '
            'due to lack of RNA reads')

    ## snps
    if args.dna_reads != '-':
        config_f= os.path.join(args.wd, 'snps', 'snps_config.yml')
        run_proc('snps', config_f, dryrun= args.dryrun, max_cores= args.cores)
    else:
        print('Single nucleotide variants will not be detected '
            'due to lack of DNA reads')
    '''
    ## denovo
    if args.dna_reads != '-':
        config_f= os.path.join(args.wd, 'denovo', 'denovo_config.yml')
        run_proc('denovo', config_f, dryrun= args.dryrun, max_cores=
                int(args.cores))
    else:
        print('The assemblies and subsequent works will not be conducted '
            'due to lack of DNA reads')

    ## phylo
    if args.dna_reads != '-':
        config_f= os.path.join(args.wd, 'phylo', 'phylo_config.yml')
        run_proc('phylo', config_f, dryrun= args.dryrun, max_cores= args.cores)
    else:
        print('The tree will not be inferred '
            'due to lack of DNA reads')
    '''
except RuntimeError as e:
    sys.exit('ERROR: {}'.format(e))
except: 
    sys.exit('Unknown problem during the analysis')
else:
    if not args.dryrun:
        try:
            ## remove redundancy
            if args.dna_reads != '-':
                config_f= os.path.join(args.wd,'cmprs_config.yml')
                run_proc('cmpr', config_f, dryrun= args.dryrun, max_cores= args.cores)
            else:
                print('The removal of redundant feature is skipped')

            ## ancestral reconstruction
            if args.dna_reads != '-' and args.rna_reads != '-':
                config_f= os.path.join(args.wd,'ar_config.yml')
                run_proc('ar', config_f, dryrun= args.dryrun, max_cores= args.cores)
            else:
                print('The ancestral reconstruction is skipped')

            # differential expression
            if args.phe_table != '-' and args.rna_reads != '-':
                config_f= os.path.join(args.wd,'de_config.yml')
                run_proc('de', config_f, dryrun= args.dryrun, max_cores= args.cores)
            else:
                print('The differential expression analysis is skipped')
        except RuntimeError as e:
            sys.exit('ERROR: {}'.format(e))
        except: 
            sys.exit('Unknown problem during the analysis')
        else:
            collect_results(args.wd)
    else:
        print('The workflow of redundancy removal, ancestral reconstruction, '
            'and differential expression analysis will be scheduled after '
            'the above processes are done')

