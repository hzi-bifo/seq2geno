#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL3

# Determine, initiate and launch the workflows based on user-defined
# arguments

import os
import sys
from tqdm import tqdm
from SGProcesses import SGProcess
import UserOptions
from CollectResults import collect_results
import LogGenerator
from Seq2GenoUtils import Warehouse
from Seq2GenoUtils import Crane
import create_config


def filter_procs(args):
    # Determine which procedures to include
    config_files = {}
    try:
        # accept config files
        if args.old_config == 'Y':
            logger.info('Skip creating config files. '
                        'Old config files will be used if found')
        else:
            config_files = create_config.main(args)
    except:
        logger.info('ERROR: fail to initiate the project')
        e = sys.exc_info()[0]
        logger.info(e)
        sys.exit()

    all_processes = []
    # >>>
    # initiate processes

    # expr
    if args.expr == 'Y':
        # ensure required data
        try:
            assert os.path.isfile(args.rna_reads)
        except:
            raise FileNotFoundError(args.rna_reads)
        all_processes.append(SGProcess(
            args.wd, 'expr',
            config_f=config_files['expr'],
            dryrun=args.dryrun,
            mem_mb=int(args.mem_mb),
            max_cores=int(args.cores)))
    else:
        logger.info('Skip counting expression levels')

    # snps
    if args.snps == 'Y':
        # ensure required data
        assert os.path.isfile(args.dna_reads)
        all_processes.append(
            SGProcess(args.wd, 'snps',
                      config_f=config_files['snps'],
                      dryrun=args.dryrun,
                      mem_mb=int(args.mem_mb),
                      max_cores=int(args.cores)))
    else:
        logger.info('Skip calling single nucleotide variants')

    # denovo
    if args.denovo == 'Y':
        # ensure required data
        assert os.path.isfile(args.dna_reads)
        all_processes.append(
            SGProcess(args.wd, 'denovo',
                      config_f=config_files['denovo'],
                      dryrun=args.dryrun,
                      mem_mb=int(args.mem_mb),
                      max_cores=int(args.cores)))
    else:
        logger.info('Skip creating de novo assemblies')

    # phylo
    if args.phylo == 'Y':
        # ensure required data
        assert os.path.isfile(args.dna_reads)
        all_processes.append(
            SGProcess(args.wd, 'phylo',
                      config_f=config_files['phylo'],
                      dryrun=args.dryrun,
                      mem_mb=int(args.mem_mb),
                      max_cores=int(args.cores)))
    else:
        logger.info('Skip inferring phylogeny')

    # ancestral reconstruction
    if args.ar == 'Y':
        # ensure required data
        assert (os.path.isfile(args.dna_reads) and
                os.path.isfile(args.rna_reads))
        all_processes.append(
            SGProcess(args.wd, 'ar',
                      config_f=config_files['ar'],
                      dryrun=args.dryrun,
                      mem_mb=int(args.mem_mb),
                      max_cores=int(args.cores)))
    else:
        logger.info('Skip ancestral reconstruction')

    # differential expression
    if args.de == 'Y':
        # ensure required data
        assert (os.path.isfile(args.phe_table) and
                os.path.isfile(args.rna_reads))
        all_processes.append(
            SGProcess(args.wd, 'de',
                      config_f=config_files['de'],
                      dryrun=args.dryrun,
                      mem_mb=int(args.mem_mb),
                      max_cores=int(args.cores)))
    else:
        logger.info('Skip differential expression analysis')

    return({'selected': all_processes, 'config_files': config_files})


def main(args):
    determined_procs = filter_procs(args)
    config_files = determined_procs['config_files']
    all_processes = determined_procs['selected']
    # >>>
    # start running processes
    try:
        processes_num = len(all_processes)+1
        pbar = tqdm(total=processes_num,
                    desc="\nseq2geno")
        for p in all_processes:
            p.run_proc()
            pbar.update(1)
    except Exception as e:
        sys.exit('ERROR: {}'.format(e))
    finally:
        if args.dryrun != 'Y':
            collect_results(args.wd, config_files)
        logger.info('Working directory {} {}'.format(
            args.wd, 'updated' if args.dryrun != 'Y' else 'unchanged'))
        pbar.update(1)

    pbar.close()


if __name__ == '__main__':
    logger = LogGenerator.make_logger()

    logger.info('Parse arguments')
    parser = UserOptions.make_parser()
    primary_args = parser.parse_args()
    # check those primary arguments
    args = UserOptions.parse_arg_yaml(primary_args.yml_f)
    args.print_args()
    # display the primary arguments only
    if primary_args.dsply_args:
        sys.exit(0)

    # determine where to run the workflows
    if primary_args.remote:
        # pack the materials and send it to the server
        new_zip_prefix = os.path.abspath(args.wd)
        new_dir = new_zip_prefix
        logger.info('Packing the materials to submit')
        Warehouse.move_data(config_f=primary_args.yml_f,
                            new_zip_prefix=new_zip_prefix,
                            new_dir=new_dir,
                            logger=logger)
        logger.info('Communicating with the server')
        sg_crane = Crane.Seq2Geno_Crane(logger=logger)
        if not os.path.isdir(args.wd):
            os.makedirs(args.wd)
        sg_crane.launch(args.wd, new_zip_prefix+'.zip')
        logger.info('DONE (remote mode)')
    else:
        # run in local machine
        main(args)
        logger.info('DONE (local mode)')
