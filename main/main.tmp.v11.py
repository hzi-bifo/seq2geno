#!/usr/bin/env python3


def main(args):
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

    ##>>>
    ## start running processes
    try:
        processes_num= len(all_processes)+1
        pbar = tqdm(total= processes_num,
            desc= "\nseq2geno")
        for p in all_processes:
            p.run_proc()
            pbar.update(1)
    except Exception as e:
        sys.exit('ERROR: {}'.format(e))
    finally:
        if not args.dryrun:
            from CollectResults import collect_results
            collect_results(args.wd, config_files)
        print('Working directory {} {}'.format(
            args.wd, 'updated' if not args.dryrun else 'unchanged'))
        pbar.update(1)

    pbar.close()

if __name__=='__main__':
    import UserOptions
    args= UserOptions.main()
    main(args)

