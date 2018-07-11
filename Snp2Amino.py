#!/usr/bin/python


"""
-------------------------------------------------------------------------------
 Name:         Snp2Amino
 Purpose:      translate SNPs to amino acids to check for substitutions

 Author:       spo12

 Created:      08.11.2013
 :q
 :q
 Changed:      28.02.2017
 Copyright:    (c) spo12 2013
-------------------------------------------------------------------------------
"""

'''    Import functions    '''
import math, argparse, sys, Bio, logging 
#import warnings
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
#from Bio import BiopythonWarning

#warnings.filterwarnings("error")

'''    Parse Options    '''

desc = """Add amino acid substitution information to a mutation table."""
epi = """\n\n"""

Parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
Parser.add_argument("-f", dest="Table", default="", metavar="<mutation table>", required=True,
                    help="path and name of the already created mutation table")
Parser.add_argument("-g", dest="GbkFile", metavar="<genbank>", required=True,
                    help="path and name of the corresponding genbank file")
Parser.add_argument("-n", dest="NonSyn", choices=['all', 'non-syn'], default="all",
                    help="specify if you want all SNPs returned or only non-synonymous ones")
Parser.add_argument("-o", dest="OutFile", metavar="<filename>", required=True,
                    help="path and name of the output file")

Args = Parser.parse_args()


'''    Global variables    '''
Table = defaultdict(list)
SnpDict = defaultdict(list)
GenDict = {}

# Set up logging so that every logging message is printed to the specified log
# file via specifying the logging level as "DEBUG".
if "." in Args.OutFile:
  name = Args.OutFile.split(".")[:-1]
  logfile = "%s.log" % "".join(name)
else:
  logfile = "%s.log" % Args.OutFile

logging.basicConfig(filename=logfile, level=logging.DEBUG, format='%(asctime)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logging.info("You have specified that you want %s SNPs." %Args.NonSyn)


'''    Read in SNP table    '''
# The table consists of the following columns: gene, pos, ref, alt and the
# names of all included isolates. The format of the gene name is:
# "locus_tag,abbreviation", but for the GBK file we only need the locus_tag.
# Every line that does not start with "gene" contains a SNP, and the only line
# that does start with "gene" is the header line.

logging.info("Reading SNP table.")
with open(Args.Table) as snp:
  for line in snp:
    if not line.startswith("gene"):
      line = line.rstrip("\n")
      info = line.split("\t")
      gene = info[0].split(",")[0] # locus
      pos = int(info[1])
      ref = info[2]
      alt = info[3]
      Table[gene].append(info)
      SnpDict[gene].append([pos, ref, alt])
    else:
      head = line


'''    Find sequences in GBK file    '''
# Extract all features from the GBK file that contains a locus_tag that appears
# as key in the SnpDict dictionary (i.e. there's a SNP in it). Find the start
# and end positions of these features to extract the sequences. Biopython can
# work directly with the start and end values.

logging.info("Reading GBK file.")
for seq_record in SeqIO.parse(Args.GbkFile, "genbank"):
  sequence = seq_record.seq
  if isinstance(sequence, Bio.Seq.UnknownSeq):
    sys.exit("There seems to be no sequence in your GenBank file!")
  for feature in seq_record.features:
    print(feature)
    '''
    ## Here should be modified!
    ## The reference .gbk  
    if not feature.type=="misc_feature" and not feature.type=="unsure": 
      if "locus_tag" in feature.qualifiers.keys():
    '''
    if (feature.type=="CDS") and ("locus_tag" in feature.qualifiers.keys()):
        locus = feature.qualifiers["locus_tag"][0]
        if locus in SnpDict.keys():
          start = feature.location.start
          end = feature.location.end
          strand = feature.location.strand
          print(strand)
          if strand == 1:
            gene_seq = seq_record.seq[start:end]
          else:
            gene_seq = seq_record.seq[start:end].reverse_complement()
          GenDict[locus] = [int(start), int(end), gene_seq, strand]
print(GenDict.keys())

'''    Exchange nucleotide and check amino acids    '''
# Loop over the genes in SnpDict. Continue if the locus_tag can also be found
# in GenDict, otherwise note that no amino acids could be exchanged.
# Loop over the SNPs for each gene and check its length. If the reference or
# the alternative is more than one nucleotide long, no amino acids can be
# exchanged.
# Next, check the strand of the gene. If the gene is on the negative strand,
# complement the SNP nucleotides. Then (this is the clever bit) substract the
# gene start position from the SNP position, just like for a plus strand SNP,
# but multiply it with -1! That way, when we replace the nucelotide at the SNP
# position via indexing, Python will start counting backwards from the end of
# the reverse-complemented sequence, which is actually the start of the
# sequence as seen from the plus strand. If the sequence is on the plus strand,
# simply substract the SNP position from the gene start, and substract 1 to
# account for Python starting to count at 0.
# We also have to figure out the sequence of the triplet containing the SNP. To
# do this, in general, we substract the start (or end) of the gene from the SNP
# position and divide the value by three. Since we don't use floats, the result
# will be rounded down. Multiplying the result with 3 again will therefore lead
# to the starting position of the triplet. If the gene is on the plus strand,
# this is all that needs to be done. For a gene on the minus strand, where we
# have to start counting triplets at the end of the gene, we're simply using the
# already negative SNP position, divide this by 3 and then multiply by 3, thanks
# to Python's rounding. The end of the triplet is then always the start +3, since
# Python uses half-open intervals.
# Once we have the in-gene SNP position, we can replace the reference with the
# alternative nucleotide by converting the string to a list of single bases and
# changing the one at the snppos index. Then we join the sequence back into one
# string, extract the triplet and translate it into an amino acid. The original
# and the mutation amino acid are then added to the SNP entry in SnpDict.

logging.info("Exchanging amino acids.")
for gene in SnpDict:
#  logging.info("%s" % GenDict[gene][2].reverse_complement())
  if gene in GenDict:
    for snp in SnpDict[gene]:
      if len(snp[1])>1 or len(snp[2])>1:
        snp.extend(("none", "none"))
      else:
        strand = GenDict[gene][3]
        logging.info("gene: %s, strand: %s" % (gene, strand))
#        logging.info("gene: %s, strand: %s, sequence: %s" % (gene, strand, GenDict[gene][2]))
        logging.info("gene fpos: %d, gene lpos: %d" % (GenDict[gene][0], GenDict[gene][1]))        
        logging.info("SNP pos: %d, reference: %s, alternative: %s" % (snp[0], snp[1], snp[2]))                
        if strand == -1:
          ref = str(Seq(snp[1]).complement())
          alt = str(Seq(snp[2]).complement())
          snppos = (snp[0] - GenDict[gene][0]) * -1
          #start = (GenDict[gene][1] + 3*((snp[0] - GenDict[gene][1]) / 3)) - GenDict[gene][0]
          start = 3 * (snppos / 3)
#          logging.info("%d - ref: %s, alt: %s, snppos: %d, start: %d" % (snp[0], ref, alt, snppos, start))
        else:
          ref = snp[1]
          alt = snp[2]
          snppos = snp[0] - GenDict[gene][0] - 1
          start = 3 * (snppos / 3)
        #try:
        end = start + 3
        gene_seq = GenDict[gene][2]
        snp_seq = list(gene_seq)
        snp_seq[snppos] = alt
        snp_seq = Seq("".join(snp_seq))
        logging.info("rel pos: %d, start: %d, end: %d, gene_seq: %s, snp_seq: %s" % (snppos, start, end, gene_seq[start:end], snp_seq[start:end]))
        if len(gene_seq[start:end]) % 3 == 0:
          orig = str(gene_seq[start:end].translate())
          muta = str(snp_seq[start:end].translate())
          snp.extend((orig, muta))
        else:
          logging.warning("Gene %s contains a partial codon. Make sure that's all right." %gene)
          snp.extend(("none", "none"))
        #except BiopythonWarning:
          #print ""gene, snp
  else:
    for snp in SnpDict[gene]:
      snp.extend(("none", "none"))


'''    Combine original table with amino acid information    '''
# Finally, the SNP table can be printed with the amino acid exchange info.
# First, add the two new column names to the header.
# Then open the output file and write the header into it, the columns separated
# by tabs. Next, check if all SNPs should be given out again, or only the non-
# synonymous ones.
# If all SNPs should be included, loop over the SnpDict dictionary, sorted by
# gene name. For each SNP in the dictionary, loop over the original mutation
# table info (Table dictionary) and find the same SNP. Add the amino acid
# exchange information to that info and write it to the output file. As soon as
# the right SNP has been found in Table break out of the last loop
# (over the Table entry).
# If only non-synonymous SNPs are wanted, begin looping over the dictionaries
# in the same way. This time, also check if the amino acid exchange information
# is not "none" for both reference and alternative, and that also reference and
# alternative are not identical. Only if this is true, add the info to the final
# line and print the SNP to the output file with it.

# Adjust header
head = head.split("\t")
head.insert(4, "ref aa")
head.insert(5, "alt aa")

logging.info("Writing output table.")
with open(Args.OutFile, "w") as out:
  out.write("\t".join(head))
  if Args.NonSyn == "all":
    for gene in sorted(SnpDict.keys()):
      for item in SnpDict[gene]:
        for line in Table[gene]:
          if item[3] == "none" or item[4] == "none" or item[3] == item[4]:
            if line[1]==str(item[0]) and line[2]==item[1] and line[3]==item[2]:
                final = line[:]
                final[4:4] = item[3:]
                out.write("\t".join(final) +"\n")
                break
  elif Args.NonSyn == "non-syn":
    for gene in sorted(SnpDict.keys()):
      for item in SnpDict[gene]:
        for line in Table[gene]:
          if item[3] != "none" and item[4] != "none" and item[3] != item[4]:
            if line[1]==str(item[0]) and line[2]==item[1] and line[3]==item[2]:
              final = line[:]
              final[4:4] = item[3:]
              out.write("\t".join(final) +"\n")
              break
