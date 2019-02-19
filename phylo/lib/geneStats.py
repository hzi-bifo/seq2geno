#!/usr/bin/env python3
from Bio import SeqIO
import os
import argparse

def run(f_l, out_d, stats_f):
    if not os.path.exists(out_d):
        os.makedirs(out_d)

    families= {}
    for line in open(f_l, 'r'):
        strain,f=line.strip().split('\t')
        print(strain)
        records= SeqIO.parse(f, 'fasta')
        for rec in records:
            family_name= str(rec.id)
            rec.id= strain
            rec.description= ''
            rec.name= ''
            if family_name in families:
                families[family_name].append(rec)
            else:
                families[family_name]=[rec]
    # write sequences
    for f in families:
        SeqIO.write(families[f], os.path.join(out_d, f+'.fa'), 'fasta')
    # statistics for furthur selection
    stats_fh= open(stats_f, 'w')
    for f in families:
        lengths= [len(rec.seq) for rec in families[f]]
        max_l= int(max(lengths))
        min_l= int(min(lengths))
        stats_fh.write('{}\t{}\t{}\n'.format(f, min_l, max_l))
    stats_fh.close()

if __name__== '__main__':
    parser= argparse.ArgumentParser(description='Compute the minimum and maximum sequences fro each family')
    parser.add_argument('--l', dest= 'f_list', required= True, help= 'list of sequence files of each strain')
    parser.add_argument('--d', dest= 'g_d', required= True, help= 'folder of output fasta files')
    parser.add_argument('--o', dest= 'stats_f', required= True, help= 'output statistics of each family')
    args= parser.parse_args()
    run(args.f_list, args.g_d, args.stats_f)

