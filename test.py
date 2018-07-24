'''
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

cores= 16
if not os.path.isdir("tmp_d/extraceted_protein_nt/"):
    os.makedirs("tmp_d/extraceted_protein_nt/")
min_samples= 2
fam_f= 'tmp/roary/gene_presence_absence.csv'

## all the files
seq_format= 'fasta'
strains= ['CH2500', 'CH2502', 'CH2522', 'F1659']
files= {strain: 'tmp/prokka/{}/{}.ffn'.format(strain, strain) for strain in strains}
print(files)

## parse the roary gpa table
fam_dict= {}
fam_df=pd.read_csv(fam_f, sep=',', quotechar='"', index_col= 0, header= 0)
sub_fam_df= fam_df[fam_df['No. isolates']>min_samples][strains]
sub_fam_df['gene']= pd.Series(sub_fam_df.index.tolist(), index= sub_fam_df.index)
sub_fam_longdf= pd.melt(sub_fam_df, id_vars=['gene'])
sub_fam_longdf= sub_fam_longdf[pd.notnull(sub_fam_longdf.value)]
#print(sub_fam_longdf.iloc[5350, :])
#print(sub_fam_longdf.value.tolist()[5350])
fam_dict= {sub_fam_longdf.value.tolist()[n]: sub_fam_longdf.gene.tolist()[n] for n in range(sub_fam_longdf.shape[0])}
#for n in range(sub_fam_longdf.shape[0]):
#    print(sub_fam_longdf.value[n])
#    print("\t"+sub_fam_longdf.gene[n])
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
    out_fasta= os.path.join("tmp_d/extraceted_protein_nt/", k+'.fa')
    out_aln= os.path.join("tmp_d/extraceted_protein_nt/", k+'.aln')
    alignments.append(out_aln)
    with open(out_fasta, 'w') as out_fh:
        out_fh.write('\n'.join(seq_dict[k]))
    aln=  MafftCommandline(quiet=True, retree= 1, thread= cores, nuc= True, globalpair= True, input= out_fasta)
    with open(out_aln, 'w') as out_fh:
        out_fh.write('\n'.join(aln()))
    list_fh.write('{}\n'.format(out_aln))
list_fh.close()
