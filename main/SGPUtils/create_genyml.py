#!/usr/bin/env python3

# Folder -> check if absent or empty
# File -> check if absent
import argparse
import yaml
import os
from LoadFile import LoadFile


def parse_usr_opts():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            description='create yaml file from Seq2Geno output for Geno2Pheno')

    parser.add_argument('-v', action='version',
                        version='v.Beta')
    parser.add_argument('--seq2geno', dest='sg',
                        help='seq2geno project folder', required=True)
    parser.add_argument('--yaml', dest='yaml',
                        help='the yaml file to be generated', required=True)

    # prediction block
    pred_args = parser.add_argument_group('predict')
    pred_args.add_argument('--opt', dest='optimize', default=['f1_macro'],
                           help='target performance metric to optimize',
                           choices=['accuracy', 'f1_pos', 'f1_macro'])
    pred_args.add_argument('--fold_n', dest='fold_n',
                           help='number of folds during validation',
                           default=10)
    pred_args.add_argument('--test_ratio', dest='test_ratio',
                           help='proportion of samples for testing',
                           default=0.1)
    pred_args.add_argument('--part', dest='part',
                           help='method to partition dataset',
                           choices=['rand'])
    pred_args.add_argument('--models', dest='models', nargs='*',
                           default=['svm'],
                           help='machine learning algorithms',
                           choices=['lsvm', 'svm', 'rf', 'lr'])
    pred_args.add_argument('--k-mer', dest='kmer', type=int,
                           help='the k-mer size for prediction',
                           default=6)
    pred_args.add_argument('--cls', dest='classes_f',
                           help='''a two-column file to specify label
                           and prediction group''')
    pred_args.add_argument('--cpu', dest='cpu', type=int, default=1,
                           help='number of cpus for parallel computation')

    args = parser.parse_args()
    return(args)


def make_genyml(args):
    blocks = dict()
    # ###
    # # metadata block
    # ###
    blocks['metadata'] = dict(
        project='sgp',
        phylogenetic_tree=os.path.abspath(os.path.join(args.sg, 'RESULTS',
                                                       'phylogeny',
                                                       'tree.nwk')),
        phenotype_table=os.path.abspath(os.path.join(args.sg, 'RESULTS',
                                                     'phenotype',
                                                     'phenotypes.mat')),
        number_of_cores=args.cpu
    )
    # ###
    # # genotype tables
    # ###
    # example:
    # dict(
    #         table= 'gpagenexp',
    #         path=
    #         "/Users/vwr33vv/Documents/project_outputs/geno2pheno/new_test/reproduce/features/features/Tobramycin_S-vs-R_features_gpa_expr.txt",
    #         preprocessing= 'none',
    #         datatype= 'text'
    #     )
    poss_tabs = {
        'kmer': dict(
            path=os.path.join(args.sg, 'RESULTS', 'assemblies'),
            preprocessing='l1',
            k_value=args.kmer),
        'gpa': dict(
            path=os.path.join(args.sg, 'RESULTS',
                              'bin_tables/gpa.mat_NONRDNT'),
            table='gpa',
            preprocessing='none',
            datatype='text'),
        'indel': dict(
            path=os.path.join(args.sg, 'RESULTS',
                              'bin_tables/indel.mat_NONRDNT'),
            table='indel',
            preprocessing='none',
            datatype='text'),
        'snp': dict(
            path=os.path.join(args.sg, 'RESULTS',
                              'bin_tables/nonsyn_SNPs_final.bin.mat_NONRDNT'),
            table='snp',
            preprocessing='none',
            datatype='text'),
        'expr': dict(
            path=os.path.join(args.sg, 'RESULTS', 'num_tables/expr.log.mat'),
            table='expr',
            preprocessing='none',
            datatype='numerical')}
    blocks['genotype_tables'] = dict(
        tables=[poss_tabs[t] for t in poss_tabs
                if (os.path.isfile(poss_tabs[t]['path'])
                    | os.path.isdir(poss_tabs[t]['path']))]
    )
    # ###
    # # prediction block
    # ###

    def parse_classes(classes_f):
        classses_dict = {}
        if args.classes_f is None:
            classses_dict = {'1': '1', '0': '0'}
        else:
            with LoadFile(classes_f) as cls_fh:
                for line in cls_fh:
                    d = line.strip().split('\t')
                    classses_dict[str(d[0])] = d[1]
        return(classses_dict)

    blocks['predictions'] = [dict(
        prediction=args.sg.strip('/').split('/')[-1],
        label_mapping=parse_classes(args.classes_f),
        optimized_for=args.optimize,
        reporting=['accuracy', 'f1_pos', 'f1_macro'],
        features=[dict(
            feature="seq2geno_feats",
            list=list(poss_tabs.keys()),
            use_validation_tuning='cv_tree')],
        classifiers={line: ("%(config)%/scikit_models/{}/"
                            "test_{}.json").format(line, line)
                     for line in args.models}
    )]
    yaml_f = args.yaml
    with LoadFile(yaml_f) as outfile:
        yaml.dump(blocks, outfile, default_flow_style=False)


if __name__ == '__main__':
    args = parse_usr_opts()
    make_genyml(args)
