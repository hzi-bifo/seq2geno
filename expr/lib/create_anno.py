#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import re
import argparse

arg_formatter = lambda prog: argparse.RawTextHelpFormatter(prog,
        max_help_position=4, width = 80)

parser = argparse.ArgumentParser(
        formatter_class= arg_formatter,
        description='create annotation table')
parser.add_argument('-r', dest= 'ref_gbk', type= str, 
    help= 'the reference genbank file')
parser.add_argument('-n', dest= 'ref_name', type= str, 
    help= 'the reference strain name')
parser.add_argument('-o', dest= 'out_f', type= str, 
    help= 'output file')

args = parser.parse_args()
#gbk_f= snakemake.input['ref_gbk']
#ref_strain= snakemake.params['ref_name']
#out_f=snakemake.output['anno_f']
gbk_f= args.ref_gbk
ref_strain= args.ref_name
out_f= args.out_f

# read the gbk
rec= SeqIO.read(open(gbk_f, 'r', encoding='windows-1252'), 'gb')
chr_len= str(len(rec.seq))
acc= rec.id
sequence_type= 'Chromosome'

with open(out_f, 'w') as out_fh:
    header_line= '@'+'|'.join([ref_strain, acc, chr_len])
    columns=['@Strain', 'Refseq_Accession', 'Replicon', 'Locus_Tag',
    'Feature_Type','Start', 'Stop', 'Strand', 'Gene_Name', 'Product_Name']
    out_fh.write(header_line+'\n')
    out_fh.write('\t'.join(columns)+'\n')
    
    target_feature= 'gene'
    for fea in [fea for fea in rec.features if fea.type == target_feature]:
        if not ('locus_tag' in fea.qualifiers):
            continue
        d=[ref_strain, acc, sequence_type, fea.qualifiers['locus_tag'][0],
            target_feature, str(fea.location.start+1), str(fea.location.end),
            '+' if fea.strand  == 1 else '-', 
            fea.qualifiers['name'][0] if 'name' in fea.qualifiers else '', 
            re.sub('\s+', '_', fea.qualifiers['function'][0]) if ('function' 
                in fea.qualifiers) else '.']
        out_fh.write('\t'.join(d)+'\n')

