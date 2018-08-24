#!/usr/bin/Rscript --vanilla 

# genes_statistics.R
#
# Perform a general statistics analysis of RNA-seq data mapped by the my_stampy_pipeline script
# Output includes an overview of the mapping results and annotation like amount of rRNA.
#
# Author: Andreas Dötsch
# last revised: Apr 2012

# parse arguments
args <- commandArgs(TRUE)
rpg_file	<- args[1]	# File containing the readcounts
anno_file	<- args[2]	# File containing the R-compatible annotation 
stats_file	<- args[3]	# File containing the Mapping statistics (output of sam_statistics.pl with -r option)
prefix		<- args[4]	# prefix for output file

# read files
rpg	<- read.table(rpg_file)
anno	<- read.table(anno_file,header=TRUE,row.names=1)
stats	<- read.table(stats_file,row.names=1)

# get global mapping stats
all_reads 	<- stats[1,]
mapped_reads 	<- stats[2,]
highQ_reads	<- stats[3,]

# count reads of gene types
all_sum  <- sum(rpg)
rRNA_sum <- sum(rpg[anno$type=="rRNA",])
sRNA_sum <- sum(rpg[anno$type=="sRNA",])
tRNA_sum <- sum(rpg[anno$type=="tRNA",])
#all_sum  <- apply(rpg,2,sum)
#rRNA_sum <- apply(rpg[anno$type=="rRNA",],2,sum)
#sRNA_sum <- apply(rpg[anno$type=="sRNA",],2,sum)
#tRNA_sum <- apply(rpg[anno$type=="tRNA",],2,sum)
mRNA_sum <- all_sum - rRNA_sum - sRNA_sum - tRNA_sum 

# calculate relative stats
mapped_perc	<- mapped_reads / all_reads * 100
highQ_perc	<- highQ_reads / all_reads * 100

rRNA_perc	<- rRNA_sum / all_sum * 100
sRNA_perc	<- sRNA_sum / all_sum * 100
tRNA_perc	<- tRNA_sum / all_sum * 100
mRNA_perc	<- mRNA_sum / all_sum * 100

# print output
out_file <- paste(prefix,".rstats",sep="")
values	<- rbind(all_reads, mapped_reads, highQ_reads, mapped_perc, highQ_perc, rRNA_sum, sRNA_sum, tRNA_sum, mRNA_sum, rRNA_perc, sRNA_perc, tRNA_perc, mRNA_perc)
write.table(values,file=out_file,quote=FALSE,sep="\t")
