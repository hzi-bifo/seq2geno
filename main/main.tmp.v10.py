#!/usr/bin/env python3

def collect_results(project_dir):
    import os
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
    if  os.path.isdir(assem_outdir):
        strains=os.listdir(assem_outdir)
        for s in strains:
            old= os.path.join(assem_outdir, s, 'contigs.fasta')
            new= os.path.join(assem_newdir, '{}.fasta'.format(s))
            if os.path.isfile(old):
                os.symlink(os.path.abspath(old), new)
    ## gpa, indel and snps tables
    bin_tab_newdir= os.path.join(results_newdir, 'bin_tables')
    os.mkdir(bin_tab_newdir)
    s_config= yaml.safe_load(open(snps_config, 'r'))
    bin_tables= [
        os.path.join(os.path.dirname(denovo_config), d_config['out_gpa_f']),
        os.path.join(os.path.dirname(denovo_config), d_config['out_indel_f']),
        os.path.join(os.path.dirname(snps_config), s_config['snps_aa_bin_mat']),
        os.path.join(os.path.dirname(snps_config), s_config['nonsyn_snps_aa_bin_mat'])
        ]
    n_rdnt_bin_tables=[f+'_NONRDNT' for f in bin_tables]
    for f in n_rdnt_bin_tables:
        new= os.path.join(bin_tab_newdir, os.path.basename(f))
        if os.path.isfile(f):
            os.symlink(os.path.abspath(f), new)
    grouping_bin_tab_newdir= os.path.join(results_newdir, 'bin_tables_groups')
    os.mkdir(grouping_bin_tab_newdir)
    grouping_bin_tables=[f+'_GROUPS' for f in bin_tables]
    for f in grouping_bin_tables:
        new= os.path.join(grouping_bin_tab_newdir, os.path.basename(f))
        if os.path.isfile(f):
            os.symlink(os.path.abspath(f), new)

    ## phylogeny
    p_config= yaml.safe_load(open(phylo_config, 'r'))
    phy_newdir= os.path.join(results_newdir, 'phylogeny')
    os.mkdir(phy_newdir)
    phylo_nwk= os.path.join(os.path.dirname(phylo_config), 
        p_config['tree_f'])
    new= os.path.join(phy_newdir, 'tree.nwk')
    if os.path.isfile(phylo_nwk):
        os.symlink(os.path.abspath(phylo_nwk), new)

    ####
    ## optional 
    ## need to check existence before collecting them
    ## expr
    e_config= yaml.safe_load(open(expr_config, 'r'))
    expr_mat= os.path.join(os.path.dirname(expr_config), 
        e_config['out_log_table'])
    if os.path.isfile(expr_mat):
        num_tab_newdir= os.path.join(results_newdir, 'num_tables')
        os.mkdir(num_tab_newdir)
        new= os.path.join(num_tab_newdir, 'expr.log.mat')
        os.symlink(os.path.abspath(expr_mat), new)

    ## ancestral reconstruction
    ar_outdir= os.path.join(os.path.dirname(expr_config), 'ancrec')
    if os.path.isdir(ar_outdir):
        ar_newdir= os.path.join(results_newdir, 'expr_ancestral_reconstructions')
        if os.path.isdir(ar_outdir):
            os.symlink(os.path.abspath(ar_outdir), ar_newdir)

    ## differential expression
    de_outdir= os.path.join(os.path.dirname(expr_config), 'dif/')
    if  os.path.isdir(de_outdir):
        de_newdir= os.path.join(results_newdir, 'differential_expression')
        if os.path.isdir(de_outdir):
            os.symlink(os.path.abspath(de_outdir), de_newdir)

    ## phenotype table
    de_config= yaml.safe_load(open(diffexpr_config, 'r'))
    phe_mat= de_config['pheno_tab']
    if os.path.isfile(phe_mat):
        phe_newdir= os.path.join(results_newdir, 'phenotype')
        os.mkdir(phe_newdir)
        new= os.path.join(phe_newdir, 'phenotypes.mat')
        if os.path.isfile(os.path.realpath(phe_mat)):
            os.symlink(os.path.abspath(os.path.realpath(phe_mat)), new)

def EditEnv(proc):
    import os
    import sys
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
    import os
    import sys
    env_dict=EditEnv(proc)
    print(proc)
    try:
        import snakemake
        import os 
        os.environ['PATH']=env_dict['PATH']
        ## unclock the snakemake working directory
        success= snakemake.snakemake(
            configfile=config_f,
            ignore_incomplete= True,
            snakefile= env_dict['SNAKEFILE'],
            workdir= os.path.dirname(config_f),
            unlock= True)
        if not success:
            raise Exception('Could not unlock target folder: {}'.format(os.path.dirname(config_f)))

        ## run the process
        success=snakemake.snakemake(
            snakefile= env_dict['SNAKEFILE'],
            restart_times= 3,
            cores= max_cores,
            configfile=config_f,
            force_incomplete= True,
            workdir= os.path.dirname(config_f),
            use_conda=True,
            conda_prefix= os.path.join(env_dict['TOOL_HOME'], 'env'),
            dryrun= dryrun,
            printshellcmds= True,
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

def main(args):
    import sys
    from tqdm import tqdm
    import os

    processes_num= 8
    pbar = tqdm(total= processes_num,
        desc= "\nseq2geno")
    try:
        ## accept config files
        if args.old_config:
            print('Skip creating config files. '
                  'Old config files will be used if found')
        else:
            import create_config
            create_config.main(args)
    except:
        print('ERROR: fail to initiate the project')
        e=sys.exc_info()[0]
        print(e)
        sys.exit()
    pbar.update(1)
    try:
        ## expr
        if args.rna_reads != '-' and args.expr:
            config_f= os.path.join(args.wd, 'expr', 'expr_config.yml')
            run_proc('expr', config_f, dryrun= args.dryrun, max_cores= args.cores)
        else:
            print('Skip counting expression levels')
        pbar.update(1)

        ## snps
        if args.dna_reads != '-' and args.snps:
            config_f= os.path.join(args.wd, 'snps', 'snps_config.yml')
            run_proc('snps', config_f, dryrun= args.dryrun, max_cores= args.cores)
        else:
            print('Skip calling single nucleotide variants')
        pbar.update(1)

        ## denovo
        if args.dna_reads != '-' and args.denovo:
            config_f= os.path.join(args.wd, 'denovo', 'denovo_config.yml')
            run_proc('denovo', config_f, dryrun= args.dryrun, max_cores=
                    int(args.cores))
        else:
            print('Skip creating de novo assemblies')
        pbar.update(1)

        ## phylo
        if args.dna_reads != '-' and args.phylo:
            config_f= os.path.join(args.wd, 'phylo', 'phylo_config.yml')
            run_proc('phylo', config_f, dryrun= args.dryrun, max_cores= args.cores)
        else:
            print('Skip inferring phylogeny')
        pbar.update(1)
    except Exception as e:
        sys.exit('ERROR: {}'.format(e))
    else:
        if not args.dryrun:
            try:
                ## ancestral reconstruction
                if args.dna_reads != '-' and args.rna_reads != '-' and args.ar:
                    config_f= os.path.join(args.wd,'ar_config.yml')
                    run_proc('ar', config_f, dryrun= args.dryrun, max_cores= args.cores)
                else:
                    print('Skip ancestral reconstruction')
                pbar.update(1)

                # differential expression
                if args.phe_table != '-' and args.rna_reads != '-' and args.de:
                    config_f= os.path.join(args.wd,'de_config.yml')
                    run_proc('de', config_f, dryrun= args.dryrun, max_cores= args.cores)
                else:
                    print('Skip differential expression analysis')
                pbar.update(1)
            except RuntimeError as e:
                sys.exit('ERROR: {}'.format(e))
            except: 
                sys.exit('ERROR: Unknown problem during the analysis')
        else:
            print('The workflow of redundancy removal, ancestral reconstruction, '
                'and differential expression analysis will be scheduled after '
                'the above processes are done')
            pbar.update(2)
    finally:
        if not args.dryrun:
            collect_results(args.wd)
        print('Working directory {} {}'.format(
            args.wd, 'updated' if not args.dryrun else 'unchanged'))
        pbar.update(1)

    pbar.close()

if __name__=='__main__':
    import UserOptions
    args= UserOptions.main()
    main(args)

