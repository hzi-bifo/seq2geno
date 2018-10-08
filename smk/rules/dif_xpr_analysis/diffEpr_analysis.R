library(DESeq2)
library(BiocParallel)
library(apeglm)

#####
#Output:
#	full DESeq2 result
#	differentially expressed genes


####
## Updates:
## Generalized for seq2geno
## Note that the input matrix of DESeq2 should have strains in COLUMNS
lib_dir= snakemake@params[['lib_dir']]
source(file.path(lib_dir, 'seq2geno.MatIO.R'))

expr_f<- snakemake@input[['expr_table']] # expression levels
target_strains<- snakemake@params[['strains']] # flat list from snakemake
phe_f<- snakemake@input[['phe_table']] # binary table of sample classes, which can include multiple columns but still follow the format of feature tables
#output_dir<- '' 
output_dir<- snakemake@output[['dif_xpr_out']]
cpu_num<- as.numeric(snakemake@params[['cores']])
alpha_cutoff<- as.numeric(snakemake@params[['alpha_cutoff']])
lfc_cutoff<- as.numeric(snakemake@params[['lfc_cutoff']])

#####
### The input matrix
### check if it in line with the feature table rules (refer to github issues and the scripts)
rpg_mat<- read_seq2geno.tab(expr_f, string_values= F)

#####
### classes
#phe_f<- '../data/pheno_table_CLSI_S-vs-R.txt'
phe_df<- read_seq2geno.tab(phe_f, string_values= T)
#colnames(phe_df)<- c('Tob', 'Cef', 'Cip', 'Mer', 'Col')
#phe_df[phe_df == 1]<- 'Resistant'
#phe_df[phe_df == 0]<- 'Sensitive'
#phe_df[is.na(phe_df)]<- 'Intermediate'
#print(phe_df[1:3, 1:3])
#print(dim(phe_df))

print(dim(rpg_mat))
print(dim(phe_df))

#####
### detect target strains

#####
### filter by row
### not necessary
#row_sums<- rowSums(rpg_mat)
#br<- c(0, 10, 100, 500, 1000, 2000, 5000, 7500, 10000, Inf)
#freq<- hist(row_sums, breaks= br, include.lowest= T, plot= F)
#ranges<- paste(head(br,-1), br[-1], sep=" - ")
#print(data.frame(range = ranges, frequency = freq$counts))
# set the cutoff= 1000
#rpg_mat<- rpg_mat[row_sums > 1000, ]

#####
### set the output directory
dir.create(output_dir, showWarnings = FALSE)
output_suffix1<- 'deseq2.tsv'
output_suffix2<- 'deseq2.DifXpr.list'

#####
### start DESeq2
for (target_col in colnames(phe_df)){
  print(target_col)
  ## detect target strains (as there may be NA phenotypes)
  target_strains<- rownames(phe_df[phe_df[,target_col] %in% c('1','0'),])
  target_strains<- target_strains[target_strains %in% rownames(rpg_mat)]
  print(target_strains)

  ## create colData
  col_df<- data.frame(sample=target_strains, pheno= phe_df[target_strains, target_col])
  rownames(col_df)<- target_strains
  dds<-DESeqDataSetFromMatrix(countData = t(rpg_mat[target_strains,]),
                              colData = col_df, 
                              design = ~pheno) 

  # determine the reference class
  dds$pheno<- factor(dds$pheno, levels= c('1', '0'))

  # start the analysis
  register(MulticoreParam(cpu_num))
  dds <- DESeq(dds, parallel= T)
  resultsNames(dds)
  res <- results(dds)
  print(dim(res))
  print(head(res))
  res_filtered <- subset(res, (padj < alpha_cutoff) & (abs(log2FoldChange)>= lfc_cutoff))
  print(dim(res_filtered))
  print(head(res_filtered))

  out_f1<- file.path(output_dir, paste0(target_col, output_suffix1, collapse= '_'))
  write.table(res, out_f1, col.names=NA, sep= '\t', quote= F)
  out_f2<- file.path(output_dir, paste0(target_col, output_suffix2, collapse= '_'))
  write(rownames(res_filtered), out_f2)

}
