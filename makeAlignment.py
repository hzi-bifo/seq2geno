'''
Create the one big concatenated alignments of each gene family
'''

from Bio import SeqIO
import textwrap
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
import os
## all the files
files= {snakemake.params['STRAINS'][n]: snakemake.input[n] for n in range(len(snakemake.params['STRAINS']))}
seq_format= 'fasta'
cores= snakemake.params['CORES']
one_big_aln_f= snakemake.output[0]
if not os.path.isdir(snakemake.params['TMP_D']):
    os.makedirs(snakemake.params['TMP_D'])

## concatenate and then split
seq_dict= {}## the values are arrays, in which the odd ones are strain names and the even ones are sequences
for s in files:
    records_dict= SeqIO.to_dict(SeqIO.parse(files[s], seq_format))
    for key in records_dict :
        if not (key in seq_dict):
            seq_dict[key]= []
        seq_dict[key].append('>'+s) # formated in fasta
        seq_dict[key].append(textwrap.fill(str(records_dict[key].seq), width= 60)) # formated in fasta

alignments= []
for k in seq_dict:
    out_fasta= os.path.join(snakemake.params['TMP_D'], k+'.fa')
    out_aln= os.path.join(snakemake.params['TMP_D'], k+'.aln')
    alignments.append(out_aln)
    with open(out_fasta, 'w') as out_fh:
        out_fh.write('\n'.join(seq_dict[k]))
    aln=  MafftCommandline(quiet=True, retree= 1, thread= cores, nuc= True, globalpair= True, input= out_fasta)
#        print(aln())
    with open(out_aln, 'w') as out_fh:
        out_fh.write('\n'.join(aln()))

one_big_aln= AlignIO.read(alignments[0], 'fasta')
one_big_aln.sort()
for f in alignments:
    aln= AlignIO.read(f, 'fasta')
    aln.sort()
    one_big_aln= one_big_aln+aln

with open(one_big_aln_f, 'w') as out_fh:
    AlignIO.write(one_big_aln, out_fh, 'fasta')
