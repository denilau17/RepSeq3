#clone composition (Fig 2B)

library(dplyr)

get.num <- function(clones){
  
  clones$cloneID <- paste0(clones$subgroup, clones$clone.num)
  
  #get cloneIDs for clones containing CD21lo cells
  cd21 <- filter(clones, cell.type=="d14.lo")
  cd21.clones <- subset(clones, clones$cloneID %in% cd21$cloneID)
  #remove Day 90 data
  cd21.clones <- filter(cd21.clones, cell.type != "d90.hi", cell.type != "d90.lo")
  
  #get the number of sequences from each cell type
  df <- cd21.clones %>% group_by(cloneID) %>% summarise(num_cells=n_distinct(cell.type))
  
  #get the number of CD21lo only cells
  CD21.only <- filter(df, num_cells==1)
  
  #get the number of CD21lo + mem clones
  CD21.2 <- filter(df, num_cells==2)
  mem.clones <- filter(clones, cell.type=="d14.hi")
  CD21.mem <- CD21.2[CD21.2$cloneID %in% mem.clones$cloneID,]
  
  #get the number of CD21lo + PB clones
  pb.clone <- filter(clones, cell.type=="d7.PB")
  CD21.pb <- CD21.2[CD21.2$cloneID %in% pb.clone$cloneID,]
  
  #get number of CD21lo + PB + Mem clones
  CD21.3 <- filter(df, num_cells==3)
  
  out <- c(nrow(CD21.only), nrow(CD21.mem), nrow(CD21.pb), nrow(CD21.3))
  return(out)
}

setwd("~/Documents/RepSeq3/clone_assignment/")
clones.007 <- read.csv("IMGT_007_85.csv", header = F)
clones.011 <- read.csv("IMGT_011_85.csv", header = F)
clones.012 <- read.csv("IMGT_012_85.csv", header = F)

header <- c("productive", "VGENE", "JGENE", "DGENE", "CDR3.len", "JUNCTION", "Sequence", "CDR3.nt",
            "CDR3.aa", "cell.type", "subject", "subgroup", "clone.num")

names(clones.007) <- header
names(clones.011) <- header
names(clones.012) <- header

get.num(clones.007)
#1487   60  118    8
get.num(clones.011)
#1991   78  217   35
get.num(clones.012)
#652  10 105 197









