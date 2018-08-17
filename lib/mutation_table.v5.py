#!/usr/bin/python
'''
Remove the loops by using pandas
Sequential coordinates are 1-based
'''

#------------------------------------------------------------------------------
# Name:         MutTab
# Purpose:      take multiple VCF files, flatcount files and the annotation to
#               create a comprehensible mutations table.
#
# Author:       spo12
#
# Created:      23.09.2013
# Last Changed: 03.06.2016
# Copyright:    (c) spo12 2013
#------------------------------------------------------------------------------

##    --        Libraries       --    ##
#Load standard, common and special libraries

import sys, re, os, time, argparse
from collections import defaultdict
import pandas as pd
from pprint import pprint

def load_flatcount(f):
    counts_str= open(f, 'r').readline().strip()
    digits= 6
    counts= [int(counts_str[i:i+digits]) for i in range(0, len(counts_str), digits)]
    counts= ['']+counts # add one at the beginning, so that the indices are 1-based
    print(len(counts))
    return(counts)
    



##    --      Parse Options     --    ##
# Ask the user for relevant files, like the "dictionary" with all files names,
# but without extensions (default is dictionary.txt) and the annotation file
# of the relevant genome (default is P. aerugiosa). Also let the user name the
# output file, if all_mutations.tab is not wanted.

desc = """Combine multiple VCF files comprehensively."""
epi = """\n\n"""

Parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
Parser.add_argument("-f", dest="DictFile", default="dictionary.txt", metavar="<dictionary>",
                    help="path and name of the dictionary file with all filenames (without extensions!)")
Parser.add_argument("-a", dest="AnnoFile", metavar="<annotation>", required=True,
                    help="path and name of the corresponding annotation file")
Parser.add_argument("-g", dest="Size", metavar="<genome size>", type=int,
                    help="(approximate) size of the reference genome, better too big than too small")
Parser.add_argument("-c", dest="Reads", metavar="<coverage/number of reads>", type=int,
                    help="minimum number of reads required at each SNP position")
Parser.add_argument("-s", dest="Score", metavar="<SNP score>", type=int,
                    help="minimum SNP score")
Parser.add_argument("-r", dest="Region", metavar="<genomic region>",
                    help="genomic region of interest (instead of whole genome), format: start-end")
Parser.add_argument("-o", dest="OutFile", default="all_mutations.tab", metavar="<filename>",
                    help="path and name of the output file")

Args = Parser.parse_args()

##    --     Global Variables   --    ##
# If no Size argument is given at the command line, the "size" variable equals
# to None, which evaluates as False (or similar) in the if clause further below.
# Initialise an empty list for isolate names and a defaultdict to store
# the mutations in.

size = Args.Size
if Args.Region:
    region = Args.Region
    first = int(region.split("-")[0])
    last = int(region.split("-")[1])
#print first, last
#print type(first), type(last) #int

names = []
mutations = defaultdict(list)

##    --     Read Annotation    --    ##
# Open the annotation file, take the first line and extract the genome size
# (usually the last information on the first line, separated by pipes "|").
# Initialise a list where each position in the genome is a 0, using the size
# extracted before.
# Skip the other header lines starting with "@" and remove trailing whitespace
# from all other lines. Split the lines at the tabs and extract locus tags
# ("gene"), gene names ("name") and gene positions. For each position in the
# genome, where the gene can be found, replace the 0 in the genome list with
# the locus tag and gene name, if applicable.

with open(Args.AnnoFile) as anno:
  if not size:
    head = anno.readline()
    finfo = head.split("|")
    try:
      size = int(finfo[-1])
    except:
      print "It seems I can't find the genome size in the annotation file!"
      print "Please enter it manually with the argument -s <number>."
  genome = [0]*size
  for line in anno:
    if not line.startswith("@"):
      line = line.rstrip()
      if line:
        info = line.split("\t")
        gene = info[3] # locus
        try:
          name = info[8] # gene name
        except:
          name = ""
        start = int(info[5])
        end = int(info[6])
        for i in range(start, end):
          if not name:
            genome[i] = gene
          else:
            genome[i] = gene +","+ name

'''
BIG CHANGE HERE
Becaue the counts are always required no matter what the parameters are, the counts are loaded first
'''
names= [f.strip() for f in open(Args.DictFile)]
names.sort()
counts_list= {}
for name in names:
    print(name)
    flc_f= 'tmp/'+name+'/bwa/variant.flatcount'
    counts= load_flatcount(flc_f)
    counts_list[name]= counts
counts_df= pd.DataFrame(data= counts_list)
print('counts_df: {}-by-{}'.format(str(counts_df.shape[0]), str(counts_df.shape[1])))

##    --      Find Mutations    --    ##
# Open the dictionary file with the relevant file names, loop over them and add
# the ".flt.vcf" extension to easily open the file. Check if the number of reads
# for each SNP should also be compared with a threshold and continue depending
# on the answer.
# If no coverage threshold is set, loop over the lines in the vcf file, skipping
# the header lines, remove trailing whitespace and extract the relevant info:
# position, reference nucleotide, alternative nucleotide, and quality. Then
# check if only a specific genomic region was selected. If not, check if a
# quality/SNP score threshold has been set. If not, save all SNPs in the
# mutations dictionary. If the SNPs should be filtered by their scores, compare
# the quality to the threshold and only add SNPs with higher scores to the
# mutations dictionary. If only a specific genomic region is of interest, first
# compare the SNP position to that region. Only if it's inside, check for the
# SNP score filter as described.
# If a coverage threshold has been set, open the ".flatcount" file as well.
# Extract the same SNP information as before, then multiply the SNP position
# by 6 (because the .flatcount file contains six digits per genome position).
# Find the right position in the .flatcount file, extract the number of reads,
# and compare that to the threshold. If the SNP position is well enough covered,
# continue with checking the genomic region and SNP score threshold as explained
# above.

#with open(Args.DictFile) as infile:
for name in names:
    filename = 'tmp/'+name+"/bwa/variant.snp-vcf"

    print(filename)
    vcf_df= pd.read_csv(filename, sep= '\t', comment= '#', 
            names= ('CHROM','POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'unknown'),
            index_col= False)
    print(vcf_df.shape)
    vcf_df['POS']= vcf_df['POS'].astype('int', copy= True)
    vcf_df['QUAL']= vcf_df['QUAL'].astype('float', copy= True)
    if Args.Region:
        vcf_df= vcf_df.loc[(vcf_df['POS'] >= first) & (vcf_df['POS'] <= last)]
    if Args.Score:
        vcf_df= vcf_df.loc[vcf_df['QUAL'] >= Args.Score]
    if Args.Reads:
        vcf_df= vcf_df[vcf_df['POS'].apply(lambda p: counts_df.loc[ p, name] > Args.Reads)]

    '''
    Big change here
    The data structure is turned into dictionary-of-dictionary
    '''
    mutations[name]= {'{}_{}_{}'.format(vcf_df['POS'][n], vcf_df['REF'][n],  vcf_df['ALT'][n]): str(vcf_df['QUAL'][n]) for n in vcf_df.index.values}

####
# the mutation dataframe
mut_df= pd.DataFrame.from_dict(data= mutations) # strains in column
print(mut_df.iloc[0:3, 0:3])
print('mut_df: {}-by-{}'.format(str(mut_df.shape[0]), str(mut_df.shape[1])))
#print(mut_df)

####
# zero coverage dataframe
print('detect zero coverage')
counts_bin_dict= {}
all_muts=mut_df.index.values.tolist()
all_muts_pos= [int(key.split("_")[0]) for key in all_muts]
counts_exp_dict= {name: {all_muts[n]: counts_df[name][all_muts_pos[n]] for n in range(len(all_muts_pos))} for name in mut_df}
counts_exp_df= pd.DataFrame.from_dict(data= counts_exp_dict)
print(counts_exp_df.iloc[0:3, 0:3])
print('counts_exp_df: {}-by-{}'.format(str(counts_exp_df.shape[0]), str(counts_exp_df.shape[1])))
#mut_df.where(mut_df.notnull(), other= 'NR', inplace= True) # if no variant and no counts, 'NR'
mut_df.where((mut_df.notnull()) | (counts_exp_df > 0), other= 'NR', inplace= True) # if no variant and no counts, 'NR'
mut_df.where((mut_df.notnull()) | (counts_exp_df == 0), other= ' ', inplace= True) # if no variant and having counts, ' '
print(mut_df.iloc[0:3, 0:3])
print('after fill-in, mut_df: {}-by-{}'.format(str(mut_df.shape[0]), str(mut_df.shape[1])))

## To use the old script
## the dataframe is transposed here
#mutations= mut_df.transpose().to_dict()
##    --       Write Output     --    ##
# Open/create the output file and write a header into it, including the names
# of the isolates. Then loop over the mutations dictionary and write each
# mutation to the file. Also include the information whether it is inside a
# gene or not by using the genome list.
#if True:# decompose the index
print('Write the table')
with open(Args.OutFile, "w") as out:
    muts= mut_df.index.values.tolist()
    muts_decomp= {mut: ['intergenic' if genome[int(mut.split('_')[0])] == 0 else genome[int(mut.split('_')[0])]] + mut.split('_') for mut in muts}
    muts_decomp_df= pd.DataFrame.from_dict(muts_decomp, orient= 'index')
    muts_decomp_df= muts_decomp_df.rename(columns={0:'gene', 1: 'pos', 2: 'ref', 3: 'alt'})
    print(muts_decomp_df.iloc[0:3, :])
    print(muts_decomp_df.shape)
    final_table= pd.concat([muts_decomp_df, mut_df], axis=1)
    print(final_table.iloc[0:3, 0:10])
    print(final_table.shape)

    final_table.to_csv(out, sep= '\t', header= True, index= False, index_label= False)
        
'''
with open(Args.OutFile, "w") as out:
  out.write("gene\tpos\tref\talt")
  for item in names:
    try:
      name = item.split("/")[-1]
    except:
      name = item
    out.write("\t"+ name)
  out.write("\n")
  for key in sorted(mutations):
    info = key.split("_")
    pos = info[0]
    ref = info[1]
    alt = info[2]
    gene = genome[int(pos)]
    if gene == 0:
      out.write("intergenic\t"+ pos +"\t"+ ref +"\t"+ alt)
    else:
      out.write(gene +"\t"+ pos +"\t"+ ref +"\t"+ alt)
    for item in sorted(mutations[key]):
      out.write("\t"+ item[1])
    out.write("\n")
'''
