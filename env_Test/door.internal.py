import os
import argparse
import subprocess

print('== in internal script ==')

parser = argparse.ArgumentParser()
functions_arg= parser.add_argument_group('functions')
functions_arg.add_argument('--env', dest='in_env', type= str, 
    default= os.getcwd(),
    help='set the working directory, where the relative paths are counted')
functions_arg.add_argument('--config', dest='in_config', type= str, 
    required= True, 
    help= 'config file for snakemake')
functions_arg.add_argument('--home', dest= 'in_h', type= str,
        required= True,
        help = 'seq2geno home')


args = parser.parse_args()
print(args)
