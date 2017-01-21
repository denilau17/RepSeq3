library(ape)

plot_tree <- function(fn){
  tree <- read.tree(fn)
  cell.type <- substr(tree$tip.label, 1,3)
  colors <- replace(cell.type, cell.type=="7PB", "dodgerblue")
  colors <- replace(colors, colors=="14l", "firebrick2")
  colors <- replace(colors, colors=="14h", "springgreen3")
  colors <- replace(colors, colors=="90h", "khaki")
  colors <- replace(colors, colors=="90l", "blueviolet")
  colors <- replace(colors, colors=="out", "black")
  tree$tip.label <- rep(NA, length(cell.type))
  plot.phylo(tree, type = "fan")
  tiplabels(pch=21, bg=colors)
}

#Figure 2E
setwd("~/Documents/RepSeq3/ml_tree/mab_tree/mab_007/")
fn_list <- list.files()
plot_tree(fn_list[4])
add.scale.bar(cex = 0.7, font = 2, col = "black")

#Figure 2F
plot_tree("~/Documents/RepSeq3/ml_tree/raxml_007/RAxML_bestTree.IGHV1-18*01_IGHJ6*02_166.fasta.phy")
add.scale.bar(cex = 0.7, font = 2, col = "black")

#Figure 2G
plot_tree("~/Documents/RepSeq3/ml_tree/raxml_012/RAxML_bestTree.IGHV1-69*01_IGHJ6*02_2213.fasta.phy")
add.scale.bar(cex = 0.7, font = 2, col = "black")


