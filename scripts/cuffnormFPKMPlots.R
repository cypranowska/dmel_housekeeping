#Title: Visualization of normalized data sets for housekeeping gene identification
#Authors: Caroline Cypranowska, Shivali Baveja
#Date: 06/21/2019

library(ggplot2)
library(reshape2)
library(dplyr)
library(gtools)

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

## join the fpkm and attr tables
# commenting out because we don't actually use this
# fpkm_attr <- merge(x = tab.use, y = gene.tab, by = "tracking_id", all.x = TRUE)

## plot the sample fpkms against each other(scatterplot)
# commenting out the line below because we don't need to create a second object to make the combos
# sample_list <- sample.key %>% pull(sample_name)
combos <- combinations(16,2,sample.key$sample_name)
# changed this for r in 2:nrow(combos), because the first column is the gene tracking_id
for (r in 2:nrow(combos)){
  # we want to plot the columns specified in combos, so we select all rows in tab.use by 
  # putting the comma first in square brackets and then choosing the column with combos[r,1] or 
  # combos[r,2]
  plot(tab.use[,combos[r,1]], tab.use[,combos[r,2]], xlab=combos[r,1], ylab=combos[r,2])
  # this will work to generate the plots but it's tricky to tell which samples are being compared in each plot
  # you can add some lines to change the labels of the x and y axis to make this more clear. I used the ggplot
  # package below, as an idea
}

# # Caroline's solution
# # convert to long format
# tab.use <- melt(tab.use)
# 
# ## generate a 'dictionary' of pairwise sample combos
# combos <- as.data.frame(combinations(16,2,sample.key$sample_name))
# # fix the levels for each column
# combos$V1 <- factor(combos$V1, levels=union(levels(combos$V1), levels(combos$V2)))
# combos$V2 <- factor(combos$V2, levels=union(levels(combos$V1), levels(combos$V2)))
# ## declare a function that makes the log10 FPKM plots
# # args: key - dataframe with pairwise sample combinations, data.use - a dataframe of log10 FPKM values for each sample in long format
# # output: .pdf file for each scatter plot
# 
# plotFPKMCombos <- function(key, data.use) {
#   for (i in 1:nrow(key)) {
#     data.use %>% filter(variable == key$V1[i] | variable == key$V2[i]) %>% dcast(tracking_id ~ variable) -> tmp
#     # remove the tracking_id column
#     tmp <- tmp[,-c(1)]
#     ggplot(tmp, aes(x = tmp[,1], y = tmp[,2])) + geom_point() + labs(x = paste("log10(FPKM+1)", colnames(tmp)[1]), y = paste("log10(FPKM+1)", colnames(tmp)[2]))
#     #ggsave(filename = paste0("log10FPKM_",colnames(tmp)[1],"_",colnames(tmp)[2],".pdf"))
#   }
# }
# 
# plotFPKMCombos(combos, tab.use)
