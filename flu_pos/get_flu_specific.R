library(dplyr)
library(tidyr)

#try to identify clones containing verfied flu positive antibodies

mab_007 <- read.csv("~/Documents/RepSeq2/flu+mAb_seqs/PB_007.csv")
mab_011 <- read.csv("~/Documents/RepSeq2/flu+mAb_seqs/PB_011.csv")
mab_012 <- read.csv("~/Documents/RepSeq2/flu+mAb_seqs/PB_012_2.csv")

setwd("~/Documents/RepSeq3/clone_assignment/")
clones_007 <- read.csv("IMGT_007_85.csv", header = F)
clones_011 <- read.csv("IMGT_011_85.csv", header = F)
clones_012 <- read.csv("IMGT_012_85.csv", header = F)

format.files <- function(df){
  names(df) <- c("productive", "VGENE", "JGENE", "DGENE", "CDR3.len", "JUNCTION", "Sequence", "CDR3.nt",
                 "CDR3.aa", "cell.type", "subject", "subgroup", "clone.num")
  df$cloneID <- paste0(df$subgroup, df$clone.num)
  df$seqID <- paste0(df$cloneID, row.names(df))
  return(df)
}

clones_007 <- format.files(clones_007)
clones_011 <- format.files(clones_011)
clones_012 <- format.files(clones_012)

get_mab_clones <- function(clones, mab){
  clones_with_mab <- clones[clones$CDR3.nt %in% mab$CDR3.IMGT,]
  all_mab_clones <- clones[clones$cloneID %in% clones_with_mab$cloneID,]
  return(all_mab_clones)
}

mab_clones_007 <- get_mab_clones(clones_007, mab_007)
mab_clones_011 <- get_mab_clones(clones_011, mab_011)
mab_clones_012 <- get_mab_clones(clones_012, mab_012)

#how many cell types are in these clones?
count.cells <- function(mab_clones){
  out <- mab_clones %>% group_by(cloneID) %>% 
    summarise(cell.types = n_distinct(cell.type), num.seqs=n())
  return(out)
} 

cells_007 <- count.cells(mab_clones_007)
cells_011 <- count.cells(mab_clones_011)
cells_012 <- count.cells(mab_clones_012)

#remove very small clones
remove_small_clones <- function(mab_clones){
  no.small.clones <- mab_clones %>% group_by(cloneID) %>% tally() %>% filter(n>=10)
  out <- mab_clones[mab_clones$cloneID %in% no.small.clones$cloneID,]
  return(out)
}

mab_clones_007 <- remove_small_clones(mab_clones_007)
mab_clones_011 <- remove_small_clones(mab_clones_011)
mab_clones_012 <- remove_small_clones(mab_clones_012)

filter_mab_clones <- function(mab){
  mab_counts <- mab %>% group_by(cloneID, cell.type) %>% summarise(count = n())
  mab_counts_wide <- spread(mab_counts, cell.type, count)
  mab_counts_wide$ratio <- mab_counts_wide$d14.lo/mab_counts_wide$d7.PB
  out <- mab_counts_wide %>% filter(d14.lo >2 & (ratio > .05 | d14.lo >= 10))
  return(out)
}

filter_mab_clones_007 <- filter_mab_clones(mab_clones_007)
filter_mab_clones_011 <- filter_mab_clones(mab_clones_011)
filter_mab_clones_012 <- filter_mab_clones(mab_clones_012)

out_007 <- mab_clones_007[mab_clones_007$cloneID %in% filter_mab_clones_007$cloneID,]
out_011 <- mab_clones_011[mab_clones_011$cloneID %in% filter_mab_clones_011$cloneID,]
out_012 <- mab_clones_012[mab_clones_012$cloneID %in% filter_mab_clones_012$cloneID,]


saveRDS(out_007, file="~/Documents/RepSeq3/flu_specific/mab_clones_007.RDS")
saveRDS(out_011, file="~/Documents/RepSeq3/flu_specific/mab_clones_011.RDS")
saveRDS(out_012, file="~/Documents/RepSeq3/flu_specific/mab_clones_012.RDS")
