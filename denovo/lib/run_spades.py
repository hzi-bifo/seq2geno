#!/usr/bin/env python3
import subprocess
import sys
out_dir= snakemake.params['spades_outdir']
cmd= []
mem_lim= int(snakemake.resources['mem_mb']/1000)
try:
    if len(snakemake.input['READS']) == 2:
        cmd= [snakemake.params['SPADES_BIN'], 
            snakemake.params['SPADES_OPT'], 
             '--threads', str(snakemake.threads), 
             '--memory', str(mem_lim),
             '-o', out_dir, 
             '-1', snakemake.input['READS'][0], 
             '-2', snakemake.input['READS'][1]]
    elif len(snakemake.input['READS']) == 1:
        cmd= [snakemake.params['SPADES_BIN'], 
            snakemake.params['SPADES_OPT'],
            '--threads', str(snakemake.threads), 
            '--memory', str(mem_lim),
             '-o', out_dir, 
            '-s', snakemake.input['READS'][0]]
    else:
        raise Exception('The reads are neither single nor paired')
    subprocess.run(cmd)
except Exception as e:
    print(e)
    sys.exit()
