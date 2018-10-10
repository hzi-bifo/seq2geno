#! /usr/bin/env python

'''
Purpose:
    Interact with the user, accept the request, and convert the information
    to the input for the subsequent snakemake workflow
'''

import os
import subprocess
import sys

def create_yaml_f(args, config_f):
    import yaml
    config_fh= open(config_f, 'w')
    yaml.dump(vars(args), config_fh, default_flow_style= False)
    config_fh.close()

if __name__== '__main__':
    seq2geno_home= os.path.abspath(
            os.path.join(os.path.dirname(os.path.realpath(__file__)),
                os.pardir))
    seq2geno_smk_dir = os.path.join(seq2geno_home, 'smk')
    seq2geno_lib_dir = os.path.join(seq2geno_home, 'lib')
    seq2geno_env_dir = os.path.join(seq2geno_home, 'env')

    # parse user's arguments
    import UserOptions
    args= UserOptions.main()
    setattr(args, 'seq2geno_smk_dir', seq2geno_smk_dir)
    setattr(args, 'seq2geno_lib_dir', seq2geno_lib_dir)
    setattr(args, 'seq2geno_env_dir', seq2geno_lib_dir)

    # create the config file
    config_f= 'config.yaml'
    create_yaml_f(args, config_f)

    # Determine which version to use (ori or ng)
    main_smk= os.path.join(seq2geno_smk_dir, 
            ('ng_MAIN.smk' if args.ng else 'MAIN.smk'))

    # Determine the main environment
    main_env= 'ng_seq2geno' if args.ng else 'seq2geno'
    main_cmd= [os.path.join(seq2geno_home, 'bin', 'BuildEnv'), 
            main_env, main_smk,
            config_f, args.workdir,
            'T' if args.dryrun else 'F', 
            'T' if args.notemp else 'F']

    # Create the environment, followed by the analysis snakemake workflows
    subprocess.call(main_cmd)
