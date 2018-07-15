'''
Because the gene names may include unsuitable characters (ie. spaces or
quotes), filenames should be 

The target gene families can be specified
The concatenation is moved to another part.
Requiring additional file to describe the families
'''
import re
import pandas as pd
from Bio import SeqIO
import textwrap
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
import os

cores= snakemake.params['CORES']
if not os.path.isdir(snakemake.params['TMP_D']):
    os.makedirs(snakemake.params['TMP_D'])
min_samples= snakemake.params['MIN_SAMPLES']
fam_f= snakemake.input['roary_gpa']
## all the files
seq_format= 'fasta'
files= {snakemake.params['STRAINS'][n]: snakemake.input['gene_dna_seqs'][n] 
        for n in range(len(snakemake.params['STRAINS']))}

## parse the roary gpa table
fam_dict= {}
fam_df=pd.read_csv(fam_f, sep=',', quotechar='"', index_col= 0, header= 0)
sub_fam_df= fam_df[fam_df['No. isolates']>min_samples][snakemake.params['STRAINS']]
sub_fam_df['gene']= pd.Series(sub_fam_df.index.tolist(), index= sub_fam_df.index)
sub_fam_longdf= pd.melt(sub_fam_df, id_vars=['gene'])
sub_fam_longdf= sub_fam_longdf[pd.notnull(sub_fam_longdf.value)]
print(sub_fam_longdf.shape)
fam_dict= {sub_fam_longdf.value.tolist()[n]: sub_fam_longdf.gene.tolist()[n]
        for n in range(sub_fam_longdf.shape[0])}
'''
for l in open(fam_f, 'r'):
    gene=re.search('^(.+):', l).group(1)
    orthologues= re.findall('(\S+)', re.search('(.+):((\s\S+)+)', l).group(2))
    for ortho in orthologues:
        fam_dict[ortho]= gene
'''
## sort by family 
## the values are arrays, in which the odd ones are strain names and the even ones are sequences
seq_dict= {}
for s in files:
    records_dict= SeqIO.to_dict(SeqIO.parse(files[s], seq_format))
    for seq_id in [x for x in records_dict.keys() if x in fam_dict] :
        key= fam_dict[seq_id] 
        if not (key in seq_dict):
            seq_dict[key]= []
        seq_dict[key].append('>'+s) # formated in fasta
        seq_dict[key].append(textwrap.fill(str(records_dict[seq_id].seq), 
            width= 60)) # formated in fasta

alignments= []
## target families
target_families= seq_dict.keys()
if 'genes_list' in snakemake.input:
    if os.path.exists(snakemake.input['genes_list']):
        target_families= [l.strip() for l in
            open(snakemake.input['genes_list'],'r')]

list_fh= open(snakemake.output['aln_list'], 'w')
for k in target_families:
    out_fasta= os.path.join(snakemake.params['TMP_D'], k+'.fa')
    out_aln= os.path.join(snakemake.params['TMP_D'], k+'.aln')
    alignments.append(out_aln)
    with open(out_fasta, 'w') as out_fh:
        out_fh.write('\n'.join(seq_dict[k]))
#    print(out_fasta)
#    print(out_aln)
#    out_fasta= re.sub('([^\w\._\-\/])', r'\\\\\1', out_fasta)
#    out_aln= re.sub('([^\w\._\-\/])', r'\\\\\1', out_aln)
    
    aln=  MafftCommandline(quiet=True, retree= 1, thread= cores, nuc= True,
            globalpair= True, input= '"'+out_fasta+'"')
    with open(out_aln, 'w') as out_fh:
        out_fh.write('\n'.join(aln()))
    list_fh.write('{}\n'.format(out_aln))
   
list_fh.close()
