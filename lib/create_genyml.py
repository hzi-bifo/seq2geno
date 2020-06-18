#!/usr/bin/env python3

'''
Folder -> check if absent or empty
File -> check if absent 
'''

def parse_usr_opts():
    import argparse
    arg_formatter = lambda prog: argparse.RawTextHelpFormatter(prog,
            max_help_position=4, width = 80)

    parser = argparse.ArgumentParser(
            formatter_class= arg_formatter,
            description='create yaml file from Seq2Geno output for Geno2Pheno')

    parser.add_argument('-v', action= 'version', 
        version='v.Beta')
    parser.add_argument('--seq2geno', dest= 'sg', 
            help= 'seq2geno project folder', required= True)
    parser.add_argument('--yaml', dest= 'yaml', 
            help= 'the yaml file to be generated', required= True)
    parser.add_argument('--out', dest= 'out', 
            help= 'output folder', default= 'geno2pheno.results')

    ## prediction block
    pred_args=  parser.add_argument_group('predict')
    pred_args.add_argument('--opt', dest= 'optimize', default= ['f1_macro'], 
            help= 'target performance metric to optimize', 
            choices=['accuracy', 'f1_pos', 'f1_macro'])
    pred_args.add_argument('--fold_n', dest= 'fold_n', 
            help= 'number of folds during validation', 
            default= 10)
    pred_args.add_argument('--test_ratio', dest= 'test_ratio', 
            help= 'proportion of samples for testing', 
            default= 0.1)
    pred_args.add_argument('--part', dest= 'part', 
            help= 'method to partition dataset', 
            choices= ['rand'])
    pred_args.add_argument('--models', dest= 'models', nargs= '*',  
            default= ['svm'],
            help= 'machine learning algorithms', 
            choices= ['lsvm', 'svm', 'rf', 'lr'])
    pred_args.add_argument('--k-mer', dest= 'kmer',type= int, 
            help= 'the k-mer size for prediction', 
            default= 6)
    pred_args.add_argument('--cls', dest= 'classes_f', 
            help= 'a two-column file to specify label and prediction group')
    pred_args.add_argument('--cpu', dest= 'cpu', type= int,default= 1,
            help= 'number of cpus for parallel computation')
    
    args = parser.parse_args()
    return(args)

def make_genyml(args):
    import yaml
    import os
    import sys

    blocks= dict()
    ####
    ## metadata block
    ####
    blocks['metadata']=dict(
        project='sgp',
        phylogenetic_tree= os.path.abspath(os.path.join(args.sg, 'RESULTS',
                                                        'phylogeny',
                                                        'tree.nwk')),
        phenotype_table= os.path.abspath(os.path.join(args.sg, 'RESULTS',
                                                      'phenotype',
                                                      'phenotypes.mat')),
        output_directory= args.out,
        number_of_cores=args.cpu
    )
    ####
    ## genotype tables
    ####
    '''
    example:
    dict( 
            table= 'gpagenexp',
            path=
            "/Users/vwr33vv/Documents/project_outputs/geno2pheno/new_test/reproduce/features/features/Tobramycin_S-vs-R_features_gpa_expr.txt",
            preprocessing= 'none',
            datatype= 'text'
        )
    '''
    poss_tabs= {
        'kmer': dict(
            path= os.path.join(args.sg,'RESULTS','assemblies'), 
            preprocessing= 'l1', 
            k_value= args.kmer),
        'gpa': dict(
            path= os.path.join(args.sg,'RESULTS','bin_tables/gpa.mat_NONRDNT'), 
            table='gpa', 
            preprocessing= 'none', 
            datatype= 'text'),
        'indel': dict(
            path= os.path.join(args.sg,'RESULTS','bin_tables/indel.mat_NONRDNT'), 
            table='indel', 
            preprocessing= 'none', 
            datatype= 'text'),
        'snp':dict(
            path= os.path.join(args.sg,'RESULTS','bin_tables/nonsyn_SNPs_final.bin.mat_NONRDNT'), 
            table='snp', 
            preprocessing= 'none', 
            datatype= 'text'),
        'expr':dict(
            path= os.path.join(args.sg, 'RESULTS', 'num_tables/expr.log.mat'), 
            table='expr', 
            preprocessing= 'none', 
            datatype= 'numerical')}
    blocks['genotype_tables']=dict(tables= 
        [poss_tabs[t] for t in poss_tabs if (os.path.isfile(poss_tabs[t]['path'])
        | os.path.isdir(poss_tabs[t]['path']))]
        )
    ####
    ## prediction block
    ####
    def parse_classes(classes_f):
        classses_dict= {}
        if args.classes_f is None:
            classses_dict= {'1': '1', '0': '0'}
        else:
            with open(classes_f, 'r') as cls_fh:
                for l in cls_fh:
                    d=l.strip().split('\t')
                    classses_dict[str(d[0])]= d[1]
        return(classses_dict)

    blocks['predictions']=[dict(
        prediction= args.sg.strip('/').split('/')[-1],
        label_mapping= parse_classes(args.classes_f),
        optimized_for= args.optimize,
        reporting= ['accuracy', 'f1_pos', 'f1_macro'], 
        features= [dict(
            feature= "seq2geno_feats",
            list= list(poss_tabs.keys()),
            use_validation_tuning=  'cv_tree')], 
#            validation_tuning= dict(name= 'cv_tree', 
#                                    train= dict(method='treebased', 
#                                               folds= args.fold_n), 
#                                    test= dict(method='treebased', 
#                                               folds= args.test_ratio), 
#                                    inner_cv= 10))
#            predefined_cv=dict(
#                name= 'predefined',
#                train=
#                "/Users/vwr33vv/Documents/project_outputs/geno2pheno/new_test/reproduce/cv/Tobramycin_S-vs-R_folds.txt", 
#                test=
#                "/Users/vwr33vv/Documents/project_outputs/geno2pheno/new_test/reproduce/cv/Tobramycin_S-vs-R_test.txt",
#                inner_cv= 10
#            )
        classifiers={l:  
                     "%(config)%/scikit_models/{}/test_{}.json".format(l, l) for l in args.models}
#        dict(
#            svm= dict(
#                param_config=
#                "/Users/vwr33vv/Documents/projects/techday/GenoPheno/data/configs/scikit_models/svm/svm_config.json")
#        )
    )]
    yaml_f= args.yaml
    with open(yaml_f, 'w') as outfile:
        yaml.dump(blocks, outfile, default_flow_style=False)

if __name__== '__main__':
    args= parse_usr_opts()
    make_genyml(args)

#    ####
#    ## genotype block
#    geno = ET.SubElement(root, "genotype")
#    #### for binary tables
#    bin_tab_dir= os.path.abspath(os.path.join(args.sg, 'RESULTS',
#        'bin_tables'))
#    if os.path.isdir(bin_tab_dir) & (len(os.listdir(bin_tab_dir) ) > 0):
#        bin_tables = ET.SubElement(
#            geno, "tables",
#            attrib= {
#                'path': bin_tab_dir, 
#                'normalization': "binary", 
#                'transpose': "False"})
#        setattr(bin_tables, 'text', 'binary_tables')
#    else:
#        raise OSError(
#                'GenML: {} (binary features) is absent but required'.format(bin_tab_dir))
#
#    #### for numeric features
#    num_tab_dir= os.path.abspath(os.path.join(args.sg, 'RESULTS',
#        'num_tables'))
#    if os.path.isdir(num_tab_dir) & (len(os.listdir(num_tab_dir) ) > 0):
#        con_tables = ET.SubElement(
#            geno, "tables",
#            attrib= {
#                'path': num_tab_dir, 
#                'normalization': "zu", 
#                'transpose': "False"})
#        setattr(con_tables, 'text', 'numeric_tables')
#    else:
#        raise OSError(
#                'GenML: {} (numeric features) is absent but required'.format(num_tables_dir))
#
#    #### genome seq
#    assem_dir= os.path.abspath(os.path.join(args.sg, 'RESULTS', 
#        'assemblies'))
#    if os.path.isdir(assem_dir) & (len(os.listdir(assem_dir) ) > 0):
#        seq= ET.SubElement(
#            geno, "sequence",
#            attrib={
#                'path': assem_dir, 
#                'kmer': str(args.kmer)})
#        setattr(seq, 'text', 'assemblies')
#    else:
#        raise OSError(
#                'GenML: {} (assemblies) is absent but required'.format(assem_dir))
#
#    ####
#    ## phenotype block
#    phe_f= os.path.abspath(os.path.join(args.sg, 'RESULTS', 'phenotype',
#        'phenotypes.mat'))
#    if os.path.isfile(phe_f):
#        pheno = ET.SubElement(
#            root, "phenotype",
#            attrib={'path': phe_f})
#        setattr(pheno, 'text', '\n')
#    else:
#        raise OSError(
#                'GenML: {} (phenotypes) is absent but required'.format(phe_f))
#
#    ####
#    ## phylogeny block
#    phy_f=os.path.abspath(os.path.join(args.sg, 'RESULTS', 'phylogeny',
#        'tree.nwk'))
#    if os.path.isfile(phy_f):
#        phy= ET.SubElement(
#            root, "phylogentictree",
#            attrib={'path': phy_f})
#        setattr(phy, 'text', '\n')
#    else:
#        raise OSError('GenML: {} (phylogeny) is absent but required'.format(phy_f))
#
#
#    ####
#    ## predict block
#    pred= ET.SubElement(root, "predict",
#            attrib= {'name': args.pred})
#    optimize= ET.SubElement(pred, "optimize")
#    setattr(optimize, 'text', str(args.optimize[0]))
#    validation= ET.SubElement(pred, "eval",
#            attrib= {'folds': str(args.fold_n), 'test': str(args.test_ratio)})
#    setattr(validation, 'text', str(args.part))
#
#    if os.path.isfile(args.classes_f):
#        all_cls= ET.SubElement(pred, "labels")
#        with open(args.classes_f, 'r') as cls_fh:
#            for l in cls_fh:
#                d=l.strip().split('\t')
#                cls= ET.SubElement(all_cls, "label",
#                        attrib= {'value': str(d[1])})
#                setattr(cls, 'text', str(d[0]))
#    else:
#        raise OSError(
#                'GenML: {} (prediction classes) is absent but required'.format(args.classes_f))
#
#    model= ET.SubElement(pred, "model")
#    if len(args.models) > 0:
#        for m in args.models:
#            ET.SubElement(model, m)
#    else:
#        raise ValueError('Prediction models not specified')
#
#    ####
#    ## convert to string
#    tree = ET.ElementTree(root)
#    rough_string= ET.tostring(root, 'utf-8')
#    reparsed= minidom.parseString(rough_string)
#    indented_string= reparsed.toprettyxml(indent="  ")
#    with open(genml_f, 'w') as genml_fh:
#            genml_fh.write(indented_string)


