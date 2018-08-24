#!/usr/bin/python

##-------------------------------------------------------------------------------
# Name:         art2cov
# Purpose:      calculate the coverage from *.art files.
#
# Author:       spo12
#
# Created:      29.05.2013
# Copyright:    (c) spo12 2013
#-------------------------------------------------------------------------------

import sys, re, os, time, subprocess, argparse


##    --      Parse Options      --    ##

#Initialise the parser, set options, parse and save them to a dictionary, then
#check if everything is there.

desc = """Calculate genome coverage of RNA-Seq data based on ART and STATS files."""
epi = """\n\n"""
#use = """usage: %prog [-h] [-f -o -u]"""

Parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
Parser.add_argument("-s", dest="Seq", metavar="PE or SE", choices=["PE", "SE"],
                help="either paired end or single end sequencing")
Parser.add_argument("-f", dest="DictFile", default="dictionary.txt", metavar="<dictionary>",
                help="path and name of the dictionary file with all filenames (without extensions!)")
Parser.add_argument("-o", dest="OutFile", default="coverage_cut2.txt", metavar="<filename>",
                help="path and name of the output file")
Parser.add_argument("-c", dest="Cov", default="2", metavar="<number>",
                help="the minimum number of reads for each position to be counted")

Args = Parser.parse_args()
if Args.Seq == "PE":
    PairedEnd = True
elif Args.Seq == "SE":
    PairedEnd = False
else:
    Parser.error("You didn't specify if you have single end or paired end data! \n Try -h if you need help.")

InFile = Args.DictFile
OutFile = file(Args.OutFile, "w")
CutOff = int(Args.Cov)


##    --     Paired End Data     --    ##

if PairedEnd:
    with open(InFile, "r") as infile:
        OutFile.write("art2cov.py using paired end data with a minimum coverage of %d reads\n\nfile\tpositions\tcoverage\ttot_mapped_reads\tmRNA" % CutOff)
        for line in infile:
            filename = line.rstrip()
            art = filename + ".art"
            stats = filename + ".rstats"
            genome_length = 0
            count = 0
            with open(art, "r") as ArtFile:
                for line in ArtFile:
                    if not line.startswith("#"):
                        genome_length += 1
                        line = line.rstrip()
                        reads = line.split("\t")
                        summ = int(reads[1]) + int(reads[2]) + int(reads[3]) + int(reads[4])
                        if summ > CutOff:
                            count += 1
            with open(stats, "r") as StatFile:
                for line in StatFile:
                    if line.startswith("mapped_reads"):
                        line = line.rstrip()
                        m_reads = line.split("\t")[1]
                    elif line.startswith("mRNA_sum"):
                        line = line.rstrip()
                        mRNA = line.split("\t")[1]
            OutFile.write("\n"+ filename +"\t"+ str(genome_length) +"\t"+
                        str(count) +"\t"+ m_reads +"\t"+ mRNA)


##    --     Single End Data     --    ##

if not PairedEnd:
    with open(InFile, "r") as infile:
        OutFile.write("art2cov.py using single end data with a minimum coverage of %d reads\n\nfile\tpositions\tcoverage\ttot_mapped_reads\tmRNA" % CutOff)
        for line in infile:
            filename = line.rstrip()
            art = filename + ".art"
            stats = filename + ".rstats"
            genome_length = 0
            count = 0
            with open(art, "r") as ArtFile:
                for line in ArtFile:
                    if not line.startswith("#"):
                        genome_length += 1
                        line = line.rstrip()
                        reads = line.split("\t")
                        summ = int(reads[1]) + int(reads[2])
                        if summ > CutOff:
                            count += 1
            with open(stats, "r") as StatFile:
                for line in StatFile:
                    if line.startswith("mapped_reads"):
                        line = line.rstrip()
                        m_reads = line.split("\t")[1]
                    elif line.startswith("mRNA_sum"):
                        line = line.rstrip()
                        mRNA = line.split("\t")[1]
            OutFile.write("\n"+ filename +"\t"+ str(genome_length) +"\t"+
                        str(count) +"\t"+ m_reads +"\t"+ mRNA)


##    --    Close Output File    --    ##

OutFile.close()
