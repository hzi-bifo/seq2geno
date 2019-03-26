import subprocess
import sys
out_dir= snakemake.params['spades_outdir']
cmd= []
## avoid simultaneously overloading
import random
import time
sleep_time= random.randrange(10)*60
time.sleep(sleep_time)

## max_trial replaced with --restart-times of snakemake
try:
    if len(snakemake.input['READS']) == 2:
        cmd= [snakemake.params['SPADES_BIN'], 
            snakemake.params['SPADES_OPT'], 
            '--threads',str(snakemake.threads), 
            '--memory', str(int(snakemake.resources['mem_mb']/1000)),
            '-o', out_dir, 
            '-1', snakemake.input['READS'][0], 
            '-2', snakemake.input['READS'][1]]
        subprocess.run(cmd)
    elif len(snakemake.input['READS']) == 1:
        cmd= [snakemake.params['SPADES_BIN'], 
            snakemake.params['SPADES_OPT'], 
            '--threads', str(snakemake.threads),
            '--memory', str(int(snakemake.resources['mem_mb']/1000)),
            '-o', out_dir, 
            '-s', snakemake.input['READS'][0]]
        subprocess.run(cmd)
    else:
        print('The reads are neither single nor paired')
        sys.exit()
except: 
    sys.exit('spades.py failed')
