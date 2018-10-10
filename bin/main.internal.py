#! /usr/bin/env python
import os
import subprocess
import argparse

def run_workflow(snakefile, config_f, workdir, dryrun, notemp):

    # run the workflow
    import snakemake

    snakemake.snakemake(
        snakefile=snakefile,
        configfile=config_f,
        unlock= True
    )
    snakemake.snakemake(
        snakefile=snakefile,
        configfile=config_f,
        workdir= workdir,
        dryrun= dryrun,
        printshellcmds= dryrun,
        notemp=notemp
    )



parser = argparse.ArgumentParser()
args= parser.add_argument_group('functions')
args.add_argument('--snakefile', dest='snakefile', type= str, required= True)
args.add_argument('--configfile', dest='configfile', type= str, required= True)
args.add_argument('--workdir', dest='workdir', type= str, required= True)
args.add_argument('--dryrun', dest='dryrun', type= str, required= True)
args.add_argument('--notemp', dest='notemp', type= str, required= True)

args = parser.parse_args()
args.dryrun= True if args.dryrun == 'T' else False
args.notemp= True if args.notemp == 'T' else False

import sys

run_workflow(args.snakefile, args.configfile, args.workdir, args.dryrun,
        args.notemp)
