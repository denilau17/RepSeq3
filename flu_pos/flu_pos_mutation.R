#flu positive mutations
library(dplyr)
library(tidyr)

mab_clones_007 <- readRDS("~/Documents/RepSeq3/flu_specific/mab_clones_007.RDS")
mab_clones_011 <- readRDS("~/Documents/RepSeq3/flu_specific/mab_clones_011.RDS")
mab_clones_012 <- readRDS("~/Documents/RepSeq3/flu_specific/mab_clones_012.RDS")

#get data from sqlite table
setwd("~/Documents/RepSeq2/")
library(RSQLite)
con <- dbConnect(RSQLite::SQLite(), "~/Documents/RepSeq2/IMGT_parsed2.sqlite")
alltables <- dbListTables(con)

IMGT_007 <- dbGetQuery( con,'select * from IMGT_007' )
IMGT_011 <- dbGetQuery( con,'select * from IMGT_011' )
IMGT_012 <- dbGetQuery( con,'select * from IMGT_012' )

IMGT_007 <- select(IMGT_007, Sequence, V.REGION.Nb.of.mutations)
IMGT_011 <- select(IMGT_011, Sequence, V.REGION.Nb.of.mutations)
IMGT_012 <- select(IMGT_012, Sequence, V.REGION.Nb.of.mutations)

IMGT_007 <- left_join(mab_clones_007, IMGT_007)
IMGT_011 <- left_join(mab_clones_011, IMGT_011)
IMGT_012 <- left_join(mab_clones_012, IMGT_012)

#check number of clones
length(unique(IMGT_007$cloneID))
length(unique(IMGT_011$cloneID))
length(unique(IMGT_012$cloneID))

#show distribution of mutation
IMGT_012$d90.lo <- NULL
IMGT_007$subject <- "007"
IMGT_011$subject <- "011"
IMGT_012$subject <- "012"

df <- rbind(IMGT_007, IMGT_011, IMGT_012)

library(ggplot2)
ggplot(df, aes(x=V.REGION.Nb.of.mutations, fill=subject)) + geom_density(alpha=0.2)

#show distributions by clone
ggplot(df, aes(x=V.REGION.Nb.of.mutations, fill=subject)) + geom_density(alpha=0.2) + 
  facet_wrap(~cloneID, ncol=2)

#plot distance between most and least mutated sequences for each clone
dist <- df %>% group_by(cloneID) %>% summarise(max=max(V.REGION.Nb.of.mutations, na.rm=T),
                                               min=min(V.REGION.Nb.of.mutations, na.rm=T),
                                               diff=max(V.REGION.Nb.of.mutations, na.rm=T)-
                                                 min(V.REGION.Nb.of.mutations, na.rm=T))
write.csv(dist, "~/Documents/RepSeq3/mutation_range.csv")

mean.mut <- function(IMGT){
  mean.df <- IMGT %>% group_by(cloneID, cell.type) %>% 
    summarise(mean.mut=mean(V.REGION.Nb.of.mutations, na.rm=T))
  mean.df <- spread(mean.df, cell.type, mean.mut)
  return(mean.df)
}

IMGT_007.mut <- mean.mut(IMGT_007)
IMGT_011.mut <- mean.mut(IMGT_011)
IMGT_012.mut <- mean.mut(IMGT_012)

IMGT_007.mut <- select (IMGT_007.mut, d14.lo, d7.PB, d14.hi)
IMGT_011.mut <- select (IMGT_011.mut, d14.lo, d7.PB, d14.hi)
IMGT_012.mut <- select (IMGT_012.mut, d14.lo, d7.PB, d14.hi)

IMGT.mut <- rbind(IMGT_007.mut, IMGT_011.mut, IMGT_012.mut)
write.csv(IMGT.mut, file="~/Documents/RepSeq3/Figure scripts/2C D/mean_mut.csv")


