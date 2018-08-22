#!/usr/bin/python

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
        gene = info[3]
        try:
          name = info[8]
        except:
          name = ""
        start = int(info[5])
        end = int(info[6])
        for i in range(start, end):
          if not name:
            genome[i] = gene
          else:
            genome[i] = gene +","+ name


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

with open(Args.DictFile) as infile:
  for line in infile:
    name = line.rstrip()
    names.append(name)
    filename = name +".flt.vcf"
    with open(filename) as vcf:
      if not Args.Reads:
        for line in vcf:
          if not line.startswith("#"):
            line = line.rstrip()
            info = line.split("\t")
            pos = info[1]
            ref = info[3]
            alt = info[4]
            qual = info[5]
            #print info
            if not Args.Region:
              if not Args.Score:
                key = pos +"_"+ ref +"_"+ alt
                mutations[key].append((name, qual))
              else:
                if float(qual)>= Args.Score:
                  key = pos +"_"+ ref +"_"+ alt
                  mutations[key].append((name, qual))
            else:
              if int(pos) >= first and int(pos) <= last:
                if not Args.Score:
                  key = pos +"_"+ ref +"_"+ alt
                  mutations[key].append((name, qual))
                else:
                  if float(qual)>= Args.Score:
                    key = pos +"_"+ ref +"_"+ alt
                    mutations[key].append((name, qual))
      else:
        flatcount = name +".flatcount"
        with open(flatcount) as flat:
          for line in vcf:
            if not line.startswith("#"):
              line = line.rstrip()
              info = line.split("\t")
              pos = info[1]
              ref = info[3]
              alt = info[4]
              qual = info[5]
              flt_pos = int(pos)*6
              flat.seek(flt_pos, 0)
              string = flat.read(6)
              if string:
                number = int(string)
                if number >= Args.Reads:
                  if not Args.Region:
                    if not Args.Score:
                      key = pos +"_"+ ref +"_"+ alt
                      mutations[key].append((name, qual))
                    else:
                      if float(qual)>= Args.Score:
                        key = pos +"_"+ ref +"_"+ alt
                        mutations[key].append((name, qual))
                  else:
                    if int(pos) >= first and int(pos) <= last:
                      if not Args.Score:
                        key = pos +"_"+ ref +"_"+ alt
                        mutations[key].append((name, qual))
                      else:
                        if float(qual)>= Args.Score:
                          key = pos +"_"+ ref +"_"+ alt
                          mutations[key].append((name, qual))

# Sort the list of samples/file names to speed up the next loop.
names.sort()


##    --     Find Read Counts   --    ##
# To check if a certain mutation does not appear in one isolate because it is
# not there and not because this region is not covered with reads, read the
# ".flatcount" file for each isolate, then loop over the mutations dictionary,
# take the position of each mutation and check the number of reads in the file.

for key in mutations:
  muts = []
  for x in mutations[key]:
    muts.append(x[0])
  for item in names:
    if item not in muts:
      filename = item +".flatcount"
      pos = int(key.split("_")[0])*6
      with open(filename) as flat:
        flat.seek(pos, 0)
        string = flat.read(6)
        if string:
          number = int(string)
          if number == 0:
            mutations[key].append((item, "NR"))
          else:
            mutations[key].append((item, " "))


##    --       Write Output     --    ##
# Open/create the output file and write a header into it, including the names
# of the isolates. Then loop over the mutations dictionary and write each
# mutation to the file. Also include the information whether it is inside a
# gene or not by using the genome list.

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
