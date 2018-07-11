
######## Snakemake header ########
import sys; sys.path.append("/home/thkuo/miniconda3/envs/phypal/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X\x16\x00\x00\x00tmp/CH2500/bwa.cons.faq\x06X\x16\x00\x00\x00tmp/CH2502/bwa.cons.faq\x07X\x16\x00\x00\x00tmp/CH2522/bwa.cons.faq\x08X\x15\x00\x00\x00tmp/F1659/bwa.cons.faq\te}q\n(X\x06\x00\x00\x00_namesq\x0b}q\x0cX\x1d\x00\x00\x00cons_coding_seqs_every_strainq\rK\x00K\x04\x86q\x0esh\rcsnakemake.io\nNamedlist\nq\x0f)\x81q\x10(h\x06h\x07h\x08h\te}q\x11h\x0b}q\x12sbubX\x06\x00\x00\x00outputq\x13csnakemake.io\nOutputFiles\nq\x14)\x81q\x15X\x0e\x00\x00\x00tmp/OneBig.alnq\x16a}q\x17(h\x0b}q\x18X\x0b\x00\x00\x00one_big_alnq\x19K\x00N\x86q\x1ash\x19h\x16ubX\x06\x00\x00\x00paramsq\x1bcsnakemake.io\nParams\nq\x1c)\x81q\x1d(K\x10X\x0c\x00\x00\x00tmp/familiesq\x1e]q\x1f(X\x06\x00\x00\x00CH2500q X\x06\x00\x00\x00CH2502q!X\x06\x00\x00\x00CH2522q"X\x05\x00\x00\x00F1659q#ee}q$(h\x0b}q%(X\x05\x00\x00\x00CORESq&K\x00N\x86q\'X\x05\x00\x00\x00TMP_Dq(K\x01N\x86q)X\x07\x00\x00\x00STRAINSq*K\x02N\x86q+uh&K\x10h(h\x1eh*h\x1fubX\t\x00\x00\x00wildcardsq,csnakemake.io\nWildcards\nq-)\x81q.X\x03\x00\x00\x00tmpq/a}q0(h\x0b}q1X\x05\x00\x00\x00TMP_Dq2K\x00N\x86q3sh(h/ubX\x07\x00\x00\x00threadsq4K\x01X\t\x00\x00\x00resourcesq5csnakemake.io\nResources\nq6)\x81q7(K\x01K\x01e}q8(h\x0b}q9(X\x06\x00\x00\x00_coresq:K\x00N\x86q;X\x06\x00\x00\x00_nodesq<K\x01N\x86q=uh:K\x01h<K\x01ubX\x03\x00\x00\x00logq>csnakemake.io\nLog\nq?)\x81q@}qAh\x0b}qBsbX\x06\x00\x00\x00configqC}qD(X\x05\x00\x00\x00tmp_dqEX\x03\x00\x00\x00tmpqFX\x08\x00\x00\x00result_dqGXL\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3/resultsqHX\x0b\x00\x00\x00working_dirqIXD\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3qJX\n\x00\x00\x00stampy_exeqKX\'\x00\x00\x00/home/thkuo/bin/stampy-1.0.23/stampy.pyqLX\t\x00\x00\x00raxml_exeqMX8\x00\x00\x00/home/thkuo/miniconda3/envs/phypal/bin/raxmlHPC-PTHREADSqNX\x11\x00\x00\x00art2genecount_exeqOX^\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/Paeru/bin/from_SusanneLab/RNA-seq/art2genecount.plqPX\x0b\x00\x00\x00sam2art_exeqQXX\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/Paeru/bin/from_SusanneLab/RNA-seq/sam2art.plqRX\x05\x00\x00\x00coresqSK\x10X\x07\x00\x00\x00samplesqTXU\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3/data/samples.tsvqUX\x12\x00\x00\x00reference_sequenceqVXh\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3/data/reference/RefCln_UCBPP-PA14.faqWX\x14\x00\x00\x00reference_annotationqXXi\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3/data/reference/RefCln_UCBPP-PA14.gbkqYX\x04\x00\x00\x00treeqZXV\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3/results/Paeru.nwkq[X\t\x00\x00\x00gpa_tableq\\XT\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3/results/gpa.tabq]X\n\x00\x00\x00expr_tableq^XU\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3/results/expr.tabq_X\t\x00\x00\x00assemblerq`X\x06\x00\x00\x00spadesqaX\t\x00\x00\x00annotatorqbX\x06\x00\x00\x00prokkaqcX\x0b\x00\x00\x00gene_sorterqdX\x05\x00\x00\x00roaryqeuX\x04\x00\x00\x00ruleqfX\x19\x00\x00\x00create_coding_regions_alnqgub.'); from snakemake.logging import logger; logger.printshellcmds = False
######## Original script #########
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
