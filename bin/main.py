#! /usr/bin/env python
import os


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
    config_f= 'test.yaml'
    create_yaml_f(args, config_f)
    print(args)

    # Determine which version to use (ori or ng)
    main_smk= os.path.join(seq2geno_smk_dir, 
            ('ng_MAIN.smk' if args.ng else 'MAIN.smk'))

    # run the workflow
    import snakemake
    snakemake.snakemake(
        snakefile=os.path.join(seq2geno_home, main_smk),
        configfile=config_f,
        unlock= True
    )
    snakemake.snakemake(
        snakefile=os.path.join(seq2geno_home, main_smk),
        configfile=config_f,
        workdir= args.workdir,
        dryrun= args.dryrun,
        notemp= args.notemp
        )

