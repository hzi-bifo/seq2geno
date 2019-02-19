#!/usr/bin/env python3
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys
from os import listdir
from os.path import join
import re
import argparse

def parse(file_list):
    ##read each alignment
    seq_d= {} ## the dictionary of all concatenated sequences
    files= [f.strip() for f in open(file_list, 'r')]
    for f in files:
        handle= open(f,'r')
        alignment= AlignIO.read(f, 'fasta')
        handle.close()
        for s in alignment:
            spe= re.match('^([^\|]+)', s.id).group(1)
            if spe in seq_d:
                seq_d[spe].append(s)
            else:
                seq_d[spe]= [s]

    return(seq_d)

def print_to_file(seq_d, out_file):
    out= open(out_file,'w')
    for spe in seq_d.keys():
        conc= seq_d[spe][0]# a record but not just a string of sequence
        conc.seq= conc.seq + ''.join([str(seq_d[spe][n].seq) for n in range(1, len(seq_d[spe]))]) 
        conc.id= spe
        conc.description= ''
        out.write(conc.format('fasta'))
    out.close()

if __name__=='__main__':
    parser= argparse.ArgumentParser(description= 'Concatenate alignments')
    parser.add_argument('--l', dest= 'aln_list', required= True, help= 'list of alignments to be concatenated')
    parser.add_argument('--o', dest= 'out', required= True, help= 'output')

    args= parser.parse_args()
    seq_groups= parse(args.aln_list)
    print_to_file(seq_groups, args.out)


