#Title: Visualization of normalized data sets for housekeeping gene identification
#Authors: Caroline Cypranowska, Shivali Baveja
#Date: 06/21/2019

library(ggplot2)
library(reshape2)
library(dplyr)

tab <- read.table("../data/classic_fpkm/genes.fpkm_table", header=T)
gene.tab <- read.table("../data/classic_fpkm/genes.attr_table", header=T)
sample.tab <- read.table("../data/classic_fpkm/samples.table",header=T)

sample.key <- read.csv("../data/classic_fpkm/sample.key.txt",header=T, stringsAsFactors = F)


## rename fpkm table columns based on sample type
names(tab) <- sample.key$sample_name[match(names(tab), sample.key$sample_id)]
names(tab)[1] <- "tracking_id"

## filter out poorly expressed genes
use <- rowSums(tab[,-1] > 5) >= 3
tab.use <- tab[use,]

## do a log10(FPKM+1) transformation of the data
tab.use[,-1] <- log10(tab.use[,-1]+1)
# convert to long format
tab.use <- melt(tab.use)

## generate a 'dictionary' of pairwise sample combos
combos <- expand.grid(sample.key$sample_name, sample.key$sample_name)
# remove rows where Var1 and Var2 are the same sample, duplicate combinations
combos <- combos[(combos$Var1 != combos$Var2),]
combos %>% mutate_if(is.factor, as.character) -> combos
combos <- combos[!duplicated(t(apply(combos,1,sort))), ]

## declare a function that makes the log10 FPKM plots
# args: key - dataframe with pairwise sample combinations, data.use - a dataframe of log10 FPKM values for each sample in long format
# output: .pdf file for each scatter plot

plotFPKMCombos <- function(key, data.use) {
  for (i in 1:nrow(key)) {
    data.use %>% filter(variable == key$Var1[i] | variable == key$Var2[i]) %>% dcast(tracking_id ~ variable) -> tmp
    # remove the tracking_id column
    tmp <- tmp[,-c(1)]
    ggplot(tmp, aes(x = tmp[,1], y = tmp[,2])) + geom_point() + labs(x = paste("log10(FPKM+1)", colnames(tmp)[1]), y = paste("log10(FPKM+1)", colnames(tmp)[2]))
    ggsave(filename = paste0("log10FPKM_",colnames(tmp)[1],"_",colnames(tmp)[2],".pdf"))
  }
}

plotFPKMCombos(combos, tab.use)
