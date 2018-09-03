library('phytools')
library('mnormt')
library('ape')
library('ggtree')

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

#args= commandArgs(trailingOnly= TRUE)
#tree_f<- args[1]
#data_f<- args[2]## two columns: names of taxon and observed values
#out_prefix<- args[3]

tree_f<- snakemake@input[['tree_f']]
data_f<- snakemake@input[['data_f']]
out_prefix<- snakemake@output[[1]]


tree<- read.newick(tree_f)
data<- as.data.frame(read_seq2geno.tab(data_f))# species in rows
genes<- colnames(data)

## check the tree
## tips without value
target_spe<- intersect(tree$tip.label, rownames(data))
if (length(tree$tip.label) < length(target_spe)){
  stop('The data do not include all tips of the tree')
}
# zero branch
if (any(tree$edge.length <= 0)){
  stop('The tree includes branch lengths equal or less than zero')
}

## run reconstruction
internal_nodes<-as.character((length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode))
tips<- tree$tip.label
out<- data.frame()
for (gene in genes){
  d<- data[,gene]
  names(d)<- row.names(data)
  rec<- fastAnc(tree= tree, x= d, vars = TRUE, CI=TRUE)
  rec_d<- rec$ace
  if (nrow(out) == 0){
    out<- data.frame(c(rec_d[internal_nodes], d[tips]))
  }
  else{
    out<- cbind(out, c(rec_d[internal_nodes], d[tips]))
  }
}
colnames(out)<- genes

## edge values
tr_info<- fortify(tree)
edges<- as.matrix(tr_info[tr_info$isTip, c('label', 'parent')])
edges<-rbind(edges, as.matrix(tr_info[! tr_info$isTip, c('node', 'parent')]))# column 1: nodes; column2: parents
edge_out<- out[edges[,1], ]-out[edges[,2], ]

out<- round(out, digits = 3)
edge_out<- round(edge_out, digits = 3)
write.table(out, file= out_prefix, quote= FALSE, sep= '\t')
write.table(edge_out, file= paste0(out_prefix, '.edge'), quote= FALSE, sep= '\t')
write.table(edges, file= paste0(out_prefix, '.parents'), quote= FALSE, row.names = FALSE, sep= '\t')
tree$node.label<- internal_nodes
write.tree(phy = tree, file = paste0(out_prefix, '.annotatedTree'))
