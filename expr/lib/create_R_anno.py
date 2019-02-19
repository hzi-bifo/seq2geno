#!/usr/bin/env python3
'''
author: tku16
group: BIFO
'''
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import re
import argparse

'''
Only the column "gene_name" will be used in the subsequently steps...

'''
arg_formatter = lambda prog: argparse.RawTextHelpFormatter(prog,
        max_help_position=4, width = 80)

parser = argparse.ArgumentParser(
        formatter_class= arg_formatter,
        description='create R annotation table')
parser.add_argument('-r', dest= 'ref_gbk', type= str, 
    help= 'the reference genbank file')
parser.add_argument('-o', dest= 'out_f', type= str, 
    help= 'output file')

#gbk_f= snakemake.input['ref_gbk']
#out_f=snakemake.output['R_anno_f']
args = parser.parse_args()
gbk_f= args.ref_gbk
out_f= args.out_f

# read the gbk
rec= ''
try:
    rec= SeqIO.read(open(gbk_f, 'r', encoding='windows-1252'), 'gb')
    chr_len= str(len(rec.seq))
    acc= rec.id
except IOError as ioe:
    import sys
    print('Unable to open {}'.format(gbk_f))
    sys.exit()

with open(out_f, 'w') as out_fh:
    columns=['locus', 'ID', 'gene_name', 'type', 'Start', 'End']
    out_fh.write('\t'.join(columns)+'\n')

    target_fea= [fea for fea in rec.features if ((fea.type == 'gene') or
        (not (re.search('RNA', fea.type) is None)))]
    for fea in [fea for fea in rec.features if ((fea.type == 'gene') or
        (not (re.search('RNA', fea.type) is None)))]:
        if not ('locus_tag' in fea.qualifiers):
            continue
        acc= fea.qualifiers['locus_tag'][0]
        gene_name= (fea.qualifiers['name'][0] if 'name' in fea.qualifiers else
                '-')
#        locus= ','.join([acc, gene_name])
        locus= ','.join([acc, fea.qualifiers['name'][0] ]) if ('name' in
        fea.qualifiers) else acc
        # be careful about the coordinates
        d=[locus, acc, gene_name, str(fea.type), str(fea.location.start+1),
                str(fea.location.end)]
        out_fh.write('\t'.join(d)+'\n')

