library(DESeq2)
library(BiocParallel)
library(apeglm)

####
## Updates:
## Generalized for seq2geno
## Note that the input matrix of DESeq2 should have strains in COLUMNS

read_seq2geno.tab<-function(f, na.value= NA){
    mat<- <- data.matrix(read.table(f, 
        sep= '\t', header= T, check.names= F, 
        stringsAsFactors= F, row.names= 1))
    if (any(is.na(mat))){
	quit()
    }else{
      return(mat)
    }
}


expr_f<- '' # expression levels
target_strains<- c() # flat list from snakemake
samples_f<- '' # binary table of sample classes, which can include multiple columns but still follow the format of feature tables
output_dir<- '' 
cpu_num<- 1
alpha_cutoff<- 0.05
lfc_cutoff<-1

#####
### The input matrix
### check if it in line with the feature table rules (refer to github issues and the scripts)
rpg_mat<- read_seq2geno.tab(expr_f, )
rpg_mat<- rpg_mat[target_strains, ]

#####
### classes
#samples_f<- '../data/pheno_table_CLSI_S-vs-R.txt'
samples_df<- read_Seq2geno.tab(samples_f, na.value= '')
#colnames(samples_df)<- c('Tob', 'Cef', 'Cip', 'Mer', 'Col')
#samples_df[samples_df == 1]<- 'Resistant'
#samples_df[samples_df == 0]<- 'Sensitive'
#samples_df[is.na(samples_df)]<- 'Intermediate'
#print(samples_df[1:3, 1:3])
#print(dim(samples_df))

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
### start DESeq2
for (target_col in colnames(samples_df)){

  ## create colData
  col_df<- data.frame(sample=target_samples, pheno= samples_df[target_samples, target_col])
  
  rownames(col_df)<- target_samples
  dds<-DESeqDataSetFromMatrix(countData = rpg_mat.T,
                              colData = col_df, 
                              design = ~pheno) 

  # determine the reference class
  dds$pheno<- relevel(dds$pheno, ref= '1')

  # start the analysis
  register(MulticoreParam(cpu_num))
  dds <- DESeq(dds, parallel= T)
  resultsNames(dds)
  res <- results(dds)
  res_filtered <- subset(res, (padj < alpha_cutoff) & (abs(log2FoldChange)>= lfc_cutoff))

  out_f<- paste(target_col, '_deseq2.txt', collapse= '')
  out_f<- paste(output_dir, out_f, collapse= '/')
  write.csv(as.data.frame(res_filtered), file=out_f)
}
