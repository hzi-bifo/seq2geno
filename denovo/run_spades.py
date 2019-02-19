import subprocess
import sys
out_dir= snakemake.params['spades_outdir']
cmd= []
max_trial= 3
for i in range(max_trial):
    print('target folder: {}'.format(out_dir))
    print('trial {}'.format(str(i)))
    try:
        if len(snakemake.input['READS']) == 2:
            cmd= [snakemake.params['SPADES_BIN'], 
                snakemake.params['SPADES_OPT'], '--threads',
                str(snakemake.threads), '-o', out_dir, 
                '-1', snakemake.input['READS'][0], 
                '-2', snakemake.input['READS'][1]]
            subprocess.run(cmd)
        elif len(snakemake.input['READS']) == 1:
            cmd= [snakemake.params['SPADES_BIN'], 
                snakemake.params['SPADES_OPT'], '--threads',
                str(snakemake.threads), '-o', out_dir, 
                '-s', snakemake.input['READS'][0]]
            subprocess.run(cmd)
        else:
            print('The reads are neither single nor paired')
            sys.exit()
    except: 
        continue
    else:
        break
