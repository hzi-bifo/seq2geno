
import subprocess
import sys
out_dir= snakemake.params['spades_outdir']
cmd= []
try:
    if len(snakemake.input['READS']) == 2:
        cmd= [snakemake.params['SPADES_BIN'], 
            snakemake.params['SPADES_OPT'], '--threads',
            str(snakemake.threads), '-o', out_dir, 
            '-1', snakemake.input['READS'][0], 
            '-2', snakemake.input['READS'][1]]
    elif len(snakemake.input['READS']) == 1:
        cmd= [snakemake.params['SPADES_BIN'], 
            snakemake.params['SPADES_OPT'], '--threads',
            str(snakemake.threads), '-o', out_dir, 
            '-s', snakemake.input['READS'][0]]
    else:
        raise Exception('The reads are neither single nor paired')
    subprocess.run(cmd)
except Exception as e:
    print(e)
    sys.exit()
