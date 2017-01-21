library(stringi)
library(ape)
library(dplyr)
library(tidyr)

#get clones with certain attributes
imgt.007 <- read.csv("~/Documents/RepSeq3/clone_assignment/IMGT_007_85.csv", header=F)
imgt.011 <- read.csv("~/Documents/RepSeq3/clone_assignment/IMGT_011_85.csv", header=F)
imgt.012 <- read.csv("~/Documents/RepSeq3/clone_assignment/IMGT_012_85.csv", header=F)

format.files <- function(df){
  names(df) <- c("productive", "VGENE", "JGENE", "DGENE", "CDR3.len", "JUNCTION", "Sequence", "CDR3.nt",
                              "CDR3.aa", "cell.type", "subject", "subgroup", "clone.num")
  df$cloneID <- paste0(df$subgroup, df$clone.num)
  df$seqID <- paste0(df$cloneID, row.names(df))
  return(df)
}

imgt.007 <- format.files(imgt.007)
imgt.011 <- format.files(imgt.011)
imgt.012 <- format.files(imgt.012)

filter_clones <- function(imgt){
  imgt_counts <- imgt %>% group_by(cloneID, cell.type) %>% summarise(count = n())
  imgt_counts_wide <- spread(imgt_counts, cell.type, count)
  imgt_counts_wide$ratio <- imgt_counts_wide$d14.lo/imgt_counts_wide$d7.PB
  out <- imgt_counts_wide %>% filter(!is.na(ratio) & d14.lo >2 & 
                                       d7.PB > 2 & (ratio > .05 | d14.lo >= 10))
  return(out)
}

imgt.007 <- filter_clones(imgt.007) 
imgt.011 <- filter_clones(imgt.011) 
imgt.012 <- filter_clones(imgt.012)

fn_007 <- paste("RAxML_bestTree",unique(imgt.007$cloneID), "fasta.phy", sep=".")
fn_011 <- paste("RAxML_bestTree",unique(imgt.011$cloneID), "fasta.phy", sep=".")
fn_012 <- paste("RAxML_bestTree",unique(imgt.012$cloneID), "fasta.phy", sep=".")

#calculates probability of identify of nearest neighbor in tree
get_probs <- function(fn, cell.list){
  #check if fn is in folder
  if (!(fn %in% dir(getwd(), "*.phy"))){
    print("no file")
    return(NULL)
  }
  tree <- read.tree(fn)
  x <- cophenetic.phylo(tree)
  x <- as.data.frame(x)
  diag(x) <- NA
  # for each row, get column index with min value. need to change 0 to NA
  #apply which.min function
  out <- apply(x, 1, function(y) which(y==min(y, na.rm=T)))
  
  #make vector of cell types
  cell.types <- names(x)
  cell.types <- substr(cell.types, 1,3)
  
  #are there multiple values?
  if (max(unlist(lapply(out, function(x) length(x))))>1){
    out <- as.data.frame(t(stri_list2matrix(out)))
    res.cell <- apply(out, 2, function(y, z) z[as.numeric(unlist(y))], z=cell.types)
    res.cell <- as.data.frame(res.cell)
  } else {
    res.cell <- as.data.frame(cell.types[as.numeric(out)])
  } 
  
  #which rows have a neighbor which is the same type of cell?
  pb <- rep(FALSE, length(cell.types))
  low <- rep(FALSE, length(cell.types))
  hi <- rep(FALSE, length(cell.types))
  
  for (i in 1:length(cell.types)){
    if ("7PB" %in% unlist(res.cell[i,])){
      pb[i] <- TRUE
    }else if ("14l" %in% unlist(res.cell[i,])){
      low[i] <- TRUE
    }else if ("14h" %in% unlist(res.cell[i,])){
      hi[i] <- TRUE
    }
  }
  
  df <- cbind.data.frame(cell.types, pb, low, hi)
  
  #calculate probabilities
  prob.lo.lo <- nrow(df %>% filter(cell.types=="14l", low==TRUE))/nrow(df %>% filter(cell.types=="14l"))
  prob.pb.lo <- nrow(df %>% filter(cell.types=="7PB", low==TRUE))/nrow(df %>% filter(cell.types=="7PB"))
  prob.hi.lo <- nrow(df %>% filter(cell.types=="14h", low==TRUE))/nrow(df %>% filter(cell.types=="14h"))
  prob.pb.pb <- nrow(df %>% filter(cell.types=="7PB", pb==TRUE))/nrow(df %>% filter(cell.types=="7PB"))
  prob.lo.pb <- nrow(df %>% filter(cell.types=="14l", pb==TRUE))/nrow(df %>% filter(cell.types=="14l"))
  prob.hi.pb <- nrow(df %>% filter(cell.types=="14h", pb==TRUE))/nrow(df %>% filter(cell.types=="14h"))
  prob.hi.hi <- nrow(df %>% filter(cell.types=="14h", hi==TRUE))/nrow(df %>% filter(cell.types=="14h"))
  prob.lo.hi <- nrow(df %>% filter(cell.types=="14l", hi==TRUE))/nrow(df %>% filter(cell.types=="14l"))
  prob.pb.hi <- nrow(df %>% filter(cell.types=="7PB", hi==TRUE))/nrow(df %>% filter(cell.types=="7PB"))
  
  df.prob <- c(prob.lo.lo, prob.pb.lo, prob.hi.lo, prob.pb.pb, prob.lo.pb, prob.hi.pb, prob.hi.hi,
               prob.lo.hi, prob.pb.hi)
  return(df.prob)
}

#all clones
setwd("~/Documents/RepSeq3/ml_tree/raxml_007/")
#fn_007 <- dir(getwd(), "*.phy")
neighbor.007 <- lapply(fn_007, get_probs, c("14l", "14h", "7PB"))
neighbor.007 <- do.call(rbind.data.frame, neighbor.007)
names(neighbor.007) <- c("low-low", "pb-low", "hi-low","pb-pb", "low-pb", 
                         "hi-pb", "hi-hi", "low-hi", "hi-pb")

setwd("~/Documents/RepSeq3/ml_tree/raxml_011/")
#fn_011 <- dir(getwd(), "*.phy")
neighbor.011 <- lapply(fn_011, get_probs, c("14l", "14h", "7PB"))
neighbor.011 <- do.call(rbind.data.frame, neighbor.011)
names(neighbor.011) <- c("low-low", "pb-low", "hi-low","pb-pb", "low-pb", 
                         "hi-pb", "hi-hi", "low-hi", "hi-pb")

setwd("~/Documents/RepSeq3/ml_tree/raxml_012/")
#fn_012 <- dir(getwd(), "*.phy")
neighbor.012 <- lapply(fn_012, get_probs, c("14l", "14h", "7PB"))
neighbor.012 <- do.call(rbind.data.frame, neighbor.012)
names(neighbor.012) <- c("low-low", "pb-low", "hi-low","pb-pb", "low-pb", 
                         "hi-pb", "hi-hi", "low-hi", "hi-pb")

csv <- rbind(neighbor.007, neighbor.011, neighbor.012)
write.csv(csv, "~/Documents/RepSeq3/ml_tree/neighbor.csv")

#mab verified clones

setwd("~/Documents/RepSeq3/ml_tree/mab_tree/mab_007/")
fn_007 <- list.files()
neighbor.007 <- lapply(fn_007, get_probs, c("14l", "14h", "7PB"))
neighbor.007 <- do.call(rbind.data.frame, neighbor.007)
names(neighbor.007) <- c("low-low", "pb-low", "hi-low","pb-pb", "low-pb", 
                         "hi-pb", "hi-hi", "low-hi", "hi-pb")

setwd("~/Documents/RepSeq3/ml_tree/mab_tree/mab_011/")
fn_011 <- list.files()
neighbor.011 <- lapply(fn_011, get_probs, c("14l", "14h", "7PB"))
neighbor.011 <- do.call(rbind.data.frame, neighbor.011)
names(neighbor.011) <- c("low-low", "pb-low", "hi-low","pb-pb", "low-pb", 
                         "hi-pb", "hi-hi", "low-hi", "hi-pb")

setwd("~/Documents/RepSeq3/ml_tree/mab_tree/mab_012/")
fn_012 <- list.files()
neighbor.012 <- lapply(fn_012, get_probs, c("14l", "14h", "7PB"))
neighbor.012 <- do.call(rbind.data.frame, neighbor.012)
names(neighbor.012) <- c("low-low", "pb-low", "hi-low","pb-pb", "low-pb", 
                         "hi-pb", "hi-hi", "low-hi", "hi-pb")

csv.flupos <- rbind(neighbor.007, neighbor.011, neighbor.012)
write.csv(csv.flupos, "~/Documents/RepSeq3/ml_tree/mab_neighbor.csv")
