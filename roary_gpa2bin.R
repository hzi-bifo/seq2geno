#!/usr/bin/Rscript

# convert the table of roary gpa (.csv) to a fasta of gpa
library('phytools')
library('ape')


load_data_matrix<- function(pa_f){
  ## with pheatmap
  # read the table
  pa_t<- read.table(pa_f, sep= ',', quote= '"', header= TRUE, row.names = 1, comment.char='', check.names=FALSE)
  pa_t<- t(pa_t)# make isolates listed in row
  heatmap<- pa_t[14:dim(pa_t)[1],] # the first row was taken as the row names
  heatmap[heatmap != '']<- '1'
  heatmap[heatmap == '']<- '0'
  storage.mode(heatmap)<-'character'
  return(heatmap)
}

pa_f<-snakemake@input[['roary_gpa']] 
print(pa_f)
out_f<-snakemake@output[['gpa_table']]
print(out_f)
mat<- load_data_matrix(pa_f)# isolates in rows

#print(dim(mat))
#print(mat[1:5, 1:5])
write.table(mat, file= out_f, sep='\t', quote= F)
