#!/usr/bin/env python3
#' Role: Manager 
#' Purpose: 
#' Determine, initiate and launch the workflows based on user-defined
#' arguments


def filter_procs(args):
    import sys
    import os
    from SGProcesses import SGProcess

    config_files= {}
    try:
        ## accept config files
        if args.old_config == 'Y':
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
    if args.expr == 'Y':
        #' ensure required data
        assert os.path.isfile(args.rna_reads)
        all_processes.append(SGProcess(
            args.wd, 'expr', 
            config_f= config_files['expr'], 
            dryrun= args.dryrun, 
            mem_mb= int(args.mem_mb), 
            max_cores=int(args.cores)))
    else:
        print('Skip counting expression levels')

    ## snps
    if args.snps == 'Y':
        #' ensure required data
        assert os.path.isfile(args.dna_reads)
        all_processes.append(SGProcess(args.wd,
                      'snps', config_f= config_files['snps'], 
                      dryrun= args.dryrun, 
                      mem_mb= int(args.mem_mb), 
                      max_cores=int(args.cores)))
    else:
        print('Skip calling single nucleotide variants')

    ## denovo
    if args.denovo == 'Y':
        #' ensure required data
        assert os.path.isfile(args.dna_reads)
        all_processes.append(SGProcess(args.wd,
                      'denovo', config_f= config_files['denovo'], 
                      dryrun= args.dryrun, 
                      mem_mb= int(args.mem_mb), 
                      max_cores=int(args.cores)))
    else:
        print('Skip creating de novo assemblies')

    ## phylo
    if args.phylo == 'Y':
        #' ensure required data
        assert os.path.isfile(args.dna_reads)
        all_processes.append(SGProcess(args.wd,
                      'phylo', config_f= config_files['phylo'], 
                      dryrun= args.dryrun, 
                      mem_mb= int(args.mem_mb), 
                      max_cores=int(args.cores)))
    else:
        print('Skip inferring phylogeny')

    ## ancestral reconstruction
    if args.ar == 'Y':
        #' ensure required data
        assert (os.path.isfile(args.dna_reads) and 
            os.path.isfile(args.rna_reads) and 
            args.dryrun != 'Y') 
        all_processes.append(SGProcess(args.wd,
                       'ar', config_f= config_files['ar'], 
                       dryrun= args.dryrun, 
                      mem_mb= int(args.mem_mb), 
                       max_cores=int(args.cores)))
    else:
        print('Skip ancestral reconstruction')

    ## differential expression
    if args.de == 'Y':
        #' ensure required data
        assert (os.path.isfile(args.phe_table) and 
            os.path.isfile(args.rna_reads) and 
            args.dryrun != 'Y') 
        all_processes.append(SGProcess(args.wd,
                      'de', config_f= config_files['de'], 
                       dryrun= args.dryrun, 
                      mem_mb= int(args.mem_mb), 
                       max_cores=int(args.cores)))
    else:
        print('Skip differential expression analysis')

    return({'selected': all_processes, 'config_files': config_files})

def main(args):
    from tqdm import tqdm
    import sys
    determined_procs= filter_procs(args)
    config_files= determined_procs['config_files']
    all_processes= determined_procs['selected']
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
        print('\n---\n')

    pbar.close()

if __name__=='__main__':
    import UserOptions
    args= UserOptions.main()
    main(args)

