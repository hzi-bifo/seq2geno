#!/usr/bin/env python3
from xml.dom import minidom
import xml.etree.ElementTree as ET
import argparse

def parse_usr_opts():
    arg_formatter = lambda prog: argparse.RawTextHelpFormatter(prog,
            max_help_position=4, width = 80)

    parser = argparse.ArgumentParser(
            formatter_class= arg_formatter,
            description='create genml file for Geno2Pheno')

    parser.add_argument('-v', action= 'version', 
        version='v.Beta')
    parser.add_argument('--genml', dest= 'genml', 
            help= 'genml filename', default= 'geno2pheno.gml')
    parser.add_argument('--project', dest= 'project', 
            help= 'project name', default= 'geno2pheno_project')
    parser.add_argument('--out', dest= 'out', 
            help= 'output folder', default= 'geno2pheno.results')

    ## genotype 
    data_args= parser.add_argument_group('data')
    data_args.add_argument('--bin_tables', dest= 'bt_dir',
        help='folder of binary features (e.g. gene pres/abs, snps...etc)',
        required= True)
    data_args.add_argument('--con_tables', dest= 'ct_dir',
        help='folder of continuous features (e.g. expression levels...etc)',
        required= True)
    data_args.add_argument('--assemblies', dest= 'assem',
        help='folder of assembled genome sequences', required= True)

    ## tree
    tree_args= parser.add_argument_group('tree')
    tree_args.add_argument('--tree', dest= 'tr',
        help='path of phylogenetic tree', required= True)

    ## phenotype block
    pheno_args= parser.add_argument_group('phenotype')
    pheno_args.add_argument('--pheno', dest= 'pheno',
        help='path of phenotype table', required= True)

    ## prediction block
    pred_args=  parser.add_argument_group('predict')
    pred_args.add_argument('--pred', dest= 'pred',
            help= 'name of prediction (e.g. classX_vs_classY)', 
            default= 'classX_vs_classY')
    pred_args.add_argument('--opt', dest= 'optimize', 
            help= 'target performance metric to optimize', 
            choices=['scores_f1_1'])
    pred_args.add_argument('--fold_n', dest= 'fold_n', 
            help= 'number of folds during validation', 
            default= 10)
    pred_args.add_argument('--test_perc', dest= 'test_perc', 
            help= 'proportion of samples for testing', 
            default= 0.1)
    pred_args.add_argument('--part', dest= 'part', 
            help= 'method to partition dataset', 
            choices= ['rand'])
    pred_args.add_argument('--models', dest= 'models', 
            help= 'machine learning algorithms', 
            choices= ['svm', 'rf', 'lr'], nargs= '*')
    pred_args.add_argument('--k-mer', dest= 'kmer', 
            help= 'the k-mer size for prediction', 
            default= 6)
    pred_args.add_argument('--cls', dest= 'classes_f', 
            help= 'a two-column file to specify label and prediction group')
    
    args = parser.parse_args()
    return(args)

if __name__== '__main__':
    args= parse_usr_opts()

    root = ET.Element("project",
            output=args.out,
            name=args.project)

    ####
    ## genotype block
    geno = ET.SubElement(root, "genotype")
    #### for binary tables
    bin_tables = ET.SubElement(geno, "bin_tables",
            path=args.bt_dir,
            normalization="binary",
            transpose="False")
    setattr(bin_tables, 'text', args.project)
    #### for numeric features
    con_tables = ET.SubElement(geno, "con_tables",
            attrib= {'path': args.ct_dir, 'normalization': "numeric",
                'transpose': "False"})
    setattr(con_tables, 'text', args.project)
    #### genome seq
    seq= ET.SubElement(geno, "sequence",
            attrib={'path': args.assem, 'kmer': str(args.kmer)})
    setattr(seq, 'text', args.project)

    ####
    ## phenotype block
    pheno = ET.SubElement(root, "phenotype",
            path=args.pheno)
    setattr(pheno, 'text', '\n')

    ####
    ## phylogeny block
    phy= ET.SubElement(root, "phylogentictree",
            path=args.tr)
    setattr(phy, 'text', '\n')

    ####
    ## predict block
    pred= ET.SubElement(root, "predict",
            attrib= {'name': args.pred})
    optimize= ET.SubElement(pred, "optimize")
    setattr(optimize, 'text', str(args.optimize))
    validation= ET.SubElement(pred, "eval",
            attrib= {'folds': str(args.fold_n), 'test': str(args.test_perc)})
    setattr(validation, 'text', str(args.part))

    with open(args.classes_f, 'r') as cls_fh:
        for l in cls_fh:
            d=l.strip().split('\t')
            cls= ET.SubElement(pred, "label",
                    attrib= {'value': str(d[0])})
            setattr(cls, 'text', str(d[1]))

    model= ET.SubElement(pred, "model")
    for m in args.models:
        ET.SubElement(model, m)


    ####
    ## convert to string
    tree = ET.ElementTree(root)
    rough_string= ET.tostring(root, 'utf-8')
    reparsed= minidom.parseString(rough_string)
    indented_string= reparsed.toprettyxml(indent="  ")
    genml_f= args.genml
    with open(genml_f, 'w') as genml_fh:
            genml_fh.write(indented_string)

