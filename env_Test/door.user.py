import os
import argparse
import subprocess

parser = argparse.ArgumentParser()
functions_arg= parser.add_argument_group('functions')
functions_arg.add_argument('--env', dest='env', type= str, 
    default= os.getcwd(),
    help='set the working directory, where the relative paths are counted')
functions_arg.add_argument('--config', dest='config', type= str, 
    required= True, 
    help= 'config file for snakemake')
functions_arg.add_argument('--home', dest= 'h', type= str,
        required= True,
        help = 'seq2geno home')


args = parser.parse_args()

print(args)
cmd= ['./set_env.sh', args.env, args.config, args.h]
subprocess.run(cmd)
