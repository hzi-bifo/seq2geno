library('phytools')
library('mnormt')
library('ape')
library('ggtree')

#####
### read the input options
script_dir<- Sys.getenv('TOOL_HOME')
source(paste(script_dir, 'seq2geno.MatIO.R', sep='/'))

tree_f<- snakemake@input[['tree_f']]
data_f<- snakemake@input[['data_f']]
out_dir<- snakemake@output[['output_dir']]
dir.create(out_dir, showWarnings = FALSE)

tree<- read.newick(tree_f)
data<- as.data.frame(read_seq2geno.tab(data_f))# species in rows
genes<- colnames(data)

#####
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
tree_backup<- tree
if (length(tree$tip.label) > length(target_spe)){
  all_tips<- tree$tip.label
  tips_to_drop<- all_tips[!all_tips %in% target_spe]
  tree<- drop.tip(tree, tip= tips_to_drop)
  pruned_tr_f<- file.path(out_dir, 'pruned_tree.nwk')
  write.tree(phy= tree, file= pruned_tr_f)
  print('The tree is pruned to fit the isolates in the phenotype data')
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

#####
### print results

out<- round(out, digits = 3)
edge_out<- round(edge_out, digits = 3)

out_node_f<- file.path(out_dir, paste(basename(data_f), 'recons_node.mat', sep= '.'))
write.table(out, file= out_node_f, quote= FALSE, sep= '\t')
out_edge_f<- file.path(out_dir, paste(basename(data_f), 'recons_edge.mat', sep= '.'))
write.table(edge_out, file= out_edge_f, quote= FALSE, sep= '\t')
out_parent_f<-file.path(out_dir, paste(basename(data_f), 'parent-to-child.tsv', sep= '.')) 
write.table(edges, file= out_parent_f, quote= FALSE, row.names = FALSE, sep= '\t')

out_id_f<-file.path(out_dir, paste(basename(tree_f), 'node_id', sep= '.')) 
tree$node.label<- internal_nodes
write.tree(phy = tree, file = out_id_f)