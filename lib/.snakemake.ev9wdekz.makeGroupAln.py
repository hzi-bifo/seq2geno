
######## Snakemake header ########
import sys; sys.path.append("/home/thkuo/miniconda3/envs/seq2geno/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X-\x00\x00\x00seq2geno_temp/roary/gene_presence_absence.csvq\x06X.\x00\x00\x00seq2geno_temp/CH2500/spades/prokka/de_novo.ffnq\x07X.\x00\x00\x00seq2geno_temp/CH2502/spades/prokka/de_novo.ffnq\x08X.\x00\x00\x00seq2geno_temp/CH2522/spades/prokka/de_novo.ffnq\tX-\x00\x00\x00seq2geno_temp/F1659/spades/prokka/de_novo.ffnq\nX.\x00\x00\x00seq2geno_temp/ESP088/spades/prokka/de_novo.ffnq\x0bX0\x00\x00\x00seq2geno_temp/MHH15083/spades/prokka/de_novo.ffnq\x0ce}q\r(X\x06\x00\x00\x00_namesq\x0e}q\x0f(X\t\x00\x00\x00roary_gpaq\x10K\x00N\x86q\x11X\r\x00\x00\x00gene_dna_seqsq\x12K\x01K\x07\x86q\x13uh\x10h\x06h\x12csnakemake.io\nNamedlist\nq\x14)\x81q\x15(h\x07h\x08h\th\nh\x0bh\x0ce}q\x16h\x0e}q\x17sbubX\x06\x00\x00\x00outputq\x18csnakemake.io\nOutputFiles\nq\x19)\x81q\x1aX=\x00\x00\x00seq2geno_temp/extracted_proteins_nt/__snakemake_dynamic__.alnq\x1ba}q\x1c(h\x0e}q\x1dX\x03\x00\x00\x00alnq\x1eK\x00N\x86q\x1fsh\x1eh\x1bubX\x06\x00\x00\x00paramsq csnakemake.io\nParams\nq!)\x81q"(K\x10K\x03]q#(X\x06\x00\x00\x00CH2500q$X\x06\x00\x00\x00CH2502q%X\x06\x00\x00\x00CH2522q&X\x05\x00\x00\x00F1659q\'X\x06\x00\x00\x00ESP088q(X\x08\x00\x00\x00MHH15083q)eX#\x00\x00\x00seq2geno_temp/extracted_proteins_ntq*e}q+(h\x0e}q,(X\x05\x00\x00\x00CORESq-K\x00N\x86q.X\x0b\x00\x00\x00MIN_SAMPLESq/K\x01N\x86q0X\x07\x00\x00\x00STRAINSq1K\x02N\x86q2X\x05\x00\x00\x00TMP_Dq3K\x03N\x86q4uh-K\x10h/K\x03h1h#h3h*ubX\t\x00\x00\x00wildcardsq5csnakemake.io\nWildcards\nq6)\x81q7X\x15\x00\x00\x00__snakemake_dynamic__q8a}q9(h\x0e}q:X\x03\x00\x00\x00famq;K\x00N\x86q<sX\x03\x00\x00\x00famq=h8ubX\x07\x00\x00\x00threadsq>K\x01X\t\x00\x00\x00resourcesq?csnakemake.io\nResources\nq@)\x81qA(K\x01K\x01e}qB(h\x0e}qC(X\x06\x00\x00\x00_coresqDK\x00N\x86qEX\x06\x00\x00\x00_nodesqFK\x01N\x86qGuhDK\x01hFK\x01ubX\x03\x00\x00\x00logqHcsnakemake.io\nLog\nqI)\x81qJ}qKh\x0e}qLsbX\x06\x00\x00\x00configqM}qN(X\x0b\x00\x00\x00working_dirqOXD\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v5qPX\n\x00\x00\x00stampy_exeqQX\'\x00\x00\x00/home/thkuo/bin/stampy-1.0.23/stampy.pyqRX\t\x00\x00\x00raxml_exeqSX\x16\x00\x00\x00raxmlHPC-PTHREADS-AVX2qTX\x08\x00\x00\x00softwareqU}qV(X\x06\x00\x00\x00mapperqWX\x03\x00\x00\x00bwaqXX\t\x00\x00\x00assemblerqYX\x06\x00\x00\x00spadesqZX\t\x00\x00\x00annotatorq[X\x06\x00\x00\x00prokkaq\\X\x0b\x00\x00\x00gene_sorterq]X\x05\x00\x00\x00roaryq^uX\x05\x00\x00\x00coresq_K\x10X\x07\x00\x00\x00samplesq`XU\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v5/data/samples.tsvqaX\x12\x00\x00\x00reference_sequenceqbXh\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v5/data/reference/RefCln_UCBPP-PA14.faqcX\x14\x00\x00\x00reference_annotationqdXi\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v5/data/reference/RefCln_UCBPP-PA14.gbkqeX\x04\x00\x00\x00treeqfXV\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v5/results/Paeru.nwkqgX\t\x00\x00\x00gpa_tableqhXT\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v5/results/gpa.tabqiX\x0e\x00\x00\x00syn_snps_tableqjXY\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v5/results/syn_snps.tabqkX\x11\x00\x00\x00nonsyn_snps_tableqlX\\\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v5/results/nonsyn_snps.tabqmX\n\x00\x00\x00expr_tableqnXU\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v5/results/expr.tabqoX\x0b\x00\x00\x00indel_tableqpXV\x00\x00\x00/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v5/results/annot.tabqquX\x04\x00\x00\x00ruleqrX\x10\x00\x00\x00sort_gene_familyqsub.'); from snakemake.logging import logger; logger.printshellcmds = False
######## Original script #########
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
