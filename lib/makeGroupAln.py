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
fam_df=pd.read_csv(fam_f, sep=',', quotechar='"', index_col= None, header= 0)
# replace the unsuitable characters
illegal_pattern= '[^\w_\-\.]'
fam_df['Gene']=fam_df['Gene'].str.replace(illegal_pattern, '-')
fam_df.set_index('Gene', inplace=True)
# reshape the data
## The paralogues in the roary output are separated, so we don't need to count
## them
sub_fam_df= fam_df[fam_df['No. isolates']>min_samples][snakemake.params['STRAINS']]
sub_fam_df['gene']= pd.Series(sub_fam_df.index.tolist(), index= sub_fam_df.index)
sub_fam_longdf= pd.melt(sub_fam_df, id_vars=['gene'])
sub_fam_longdf= sub_fam_longdf[pd.notnull(sub_fam_longdf.value)]
print(sub_fam_longdf.shape)
## gene name -> gene family
fam_dict= {sub_fam_longdf.value.tolist()[n]: sub_fam_longdf.gene.tolist()[n]
        for n in range(sub_fam_longdf.shape[0])}
'''
for l in open(fam_f, 'r'):
    gene=re.search('^(.+):', l).group(1)
    orthologues= re.findall('(\S+)', re.search('(.+):((\s\S+)+)', l).group(2))
    for ortho in orthologues:
        fam_dict[ortho]= gene
'''
## obtain sequence information, sort them by family, and format the seqeuences
seq_dict= {}
for strain in files:
    records_dict= SeqIO.to_dict(SeqIO.parse(files[strain], seq_format))
    # extract the seqeunce and format it
    seq_series= sub_fam_df[strain].apply(lambda x:
            '\n'.join(['>'+strain, textwrap.fill(str(records_dict[x].seq),
                width= 60), '']) if x in records_dict else
            '\n'.join(['>'+strain,'']))
    seq_dict[strain]= seq_series
seq_df= pd.DataFrame(data=seq_dict) # strains in columns
seq_df=seq_df.T # strains in rows
print(seq_df.shape)

alignments= []
## target families
#target_families= seq_dict.keys()
target_families= seq_df.columns.values.tolist()
if 'genes_list' in snakemake.input:
    if os.path.exists(snakemake.input['genes_list']):
        target_families= [l.strip() for l in
            open(snakemake.input['genes_list'],'r')]

## The list is replaced with the 'dynamic' function of snakemake
#list_fh= open(snakemake.output['aln_list'], 'w')
for k in target_families:
    out_fasta= os.path.join(snakemake.params['TMP_D'], k+'.fa')
    out_aln= os.path.join(snakemake.params['TMP_D'], k+'.aln')
    out_fasta= re.sub('([^\w_\-\/\.])', r'\\\1', out_fasta)
    out_aln= re.sub('([^\w_\-\/\.])', r'\\\1', out_aln)
    alignments.append(out_aln)
    with open(out_fasta, 'w') as out_fh:
        out_fh.write('\n'.join(seq_df[k].values.tolist()))
        #out_fh.write('\n'.join(seq_dict[k]))
#    print(out_fasta)
#    print(out_aln)
    
    aln=  MafftCommandline(quiet=True, retree= 1, thread= cores, nuc= True,
            globalpair= True, input= '"'+out_fasta+'"')
    with open(out_aln, 'w') as out_fh:
        out_fh.write('\n'.join(aln()))
#    list_fh.write('{}\n'.format(out_aln))
   
#list_fh.close()
