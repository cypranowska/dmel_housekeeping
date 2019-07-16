#Title: Cleaning Sample Data in order to prevent across the board low-expressing genes from skewing the data.
#Authors: Caroline Cypranowska, Shivali Baveja
#Date: 07/09/2019

###from cuffnormFPKMPlots.R
## import data for cleaning
tab <- read.table("../data/classic_fpkm/genes.fpkm_table", header=T)
sample.key <- read.csv("../data/classic_fpkm/sample.key.txt",header=T, stringsAsFactors = F)

## rename fpkm table columns based on sample type
names(tab) <- sample.key$sample_name[match(names(tab), sample.key$sample_id)]
names(tab)[1] <- "tracking_id"

## filter out poorly expressed genes
use <- rowSums(tab[,-1] > 5) >= 3
tab.use <- tab[use,]

## do a log10(FPKM+1) transformation of the data
tab.use[,-1] <- log10(tab.use[,-1]+1)

##running ANOVA 
for (gene in 2:nrows(tab.use)){
  fit <- aov(tab.use[gene,] ~ tracking_id)
  summary(fit, test="Pillai")
  if (fit$Pr(>F) <= 0.05) {
    tab.use <- tab.use[-c(gene),]}
}

## sort selected genes by highest average expression across tissue types
tab.averages <- data.frame(gene=tab.use[,1], avg_exp=rowMeans(tab.use[,-1]))
tab.sorted <- tab.averages[with(tab.averages,order(-avg_exp)),]

##view 100 genes with highest avg expression
head(tab.sorted, n=100)

## plot the sample fpkms against each other(scatterplot)
combos <- combinations(16,2,sample.key$sample_name)
for (r in 2:nrow(combos)){
  plot(tab.use[,combos[r,1]], tab.use[,combos[r,2]], main=cat(sprintf("%s vs %s", combos[r,1], combos[r,2])), xlab=combos[r,1],ylab=combos[r,2])
}

