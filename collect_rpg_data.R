#!/usr/bin/Rscript --vanilla 

# collect_rpg_data.R
#
# Reads all data from rpg gene expression files specified in a dictionary file and writes it to a single table, which can be used e.g. for DESeq.
#
# Author: Andreas Dötsch
# last revised: Sep 2012

# parse arguments
#args <- commandArgs(TRUE)
#dict_file	<- args[1]	# File containing the sample dictionary
#anno_file	<- args[2]	# annotations of genes
dict_file	<- snakemake@input[['dict_f']]	# File containing the sample dictionary
anno_file	<- snakemake@params[['ANNO_F']]	# annotations of genes

# read annotation
#annotation <- read.table("/data/reference_sequences/Pseudomonas_aeruginosa_PA14_R_annotation",header=TRUE,row.names=1)
#annotation <- read.table("/data/reference_sequences/Pseudomonas_aeruginosa_PAO1_R_annotation",header=TRUE,row.names=1)
#annotation <- read.table("R_annotation_Denitsas_TSS.tab",header=TRUE,row.names=1)
annotation <- read.table(anno_file,header=TRUE,row.names=1)

# set path
#path <- "/data2/RNA-seq/transcriptome_c-di-GMP/01_2013/"
#path <- "/data2/RNA-seq/"
#path <- ""

# read & parse dictionary
dict			<- read.table(dict_file,header=TRUE, stringsAsFactor= F, check.names= F)
nsamples		<- nrow(dict)			# no. of samples
all_rpg			<- matrix(0,length(annotation$gene_name),nsamples)
rownames(all_rpg)	<- rownames(annotation)
colnames(all_rpg)	<- dict$name
#conditions		<- dict$condition

#for(i in 1:nsamples){
#	rpg_file	<- as.character(dict[i,3])	# path to rpg file
#	stats_file	<- dict[i,4]	# path to rstats file	=> TODO: implement stats for quality control

#print(i)
for (i in 1:nrow(dict)){
	# parse rpg file
#	rpg		<- as.matrix(read.table(paste(path,rpg_file,sep="")))
  print(dict[i, ]$path)
	rpg		<- data.matrix(read.table(dict[i, ]$path))
  print(length(rpg))
	if(ncol(rpg) > 1){	
		# two-strand mode (default setting of art2genecount.pl)
		rpg	<- as.numeric(rpg[,2])
	} # otherwise blank mode is assumed
	all_rpg[,i]	<- rpg
}
write.table(all_rpg,paste(dict_file,".rpg",sep=""),sep="\t", col.names=NA, quote= F)
