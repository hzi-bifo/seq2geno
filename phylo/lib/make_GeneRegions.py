#!/usr/bin/env python2
import sys
import re
import argparse

def select(f, target_feature, comment_pattern):
    lines=[l.strip().split('\t') for l in open(f, 'r') if (not re.search(comment_pattern,l))&(len(l.split('\t')) == 9)]
    regions= ['{}:{}-{}'.format(l[0], l[3], l[4]) for l in lines if l[2] == target_feature]

    if len(regions) == 0:
	sys.exit('No target lines found')
    else:
	print('\n'.join(regions))

if __name__=='__main__':
    
    parser= argparse.ArgumentParser(description='Create gene regions list')
    parser.add_argument('--g', dest= 'gff', required= True, help= 'gff file')
    parser.add_argument('--f', dest= 'feature', required= True, help= 'target features (column 3 of gff file)')
    parser.add_argument('--c', dest= 'comment', required= False, default= '^#', help= 'comment pattern')
    args= parser.parse_args()

    f= args.gff
    target_feature= args.feature
    comment_pattern= args.comment
    select(f, target_feature, comment_pattern)
