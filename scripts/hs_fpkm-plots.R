

tab <- read.table("../data/classic_fpkm/genes.fpkm_table", header=T)
gene.tab <- read.table("../data/classic_fpkm/genes.attr_table", header=T)
sample.tab <- read.table("../data/classic_fpkm/samples.table",header=T)

sample.key <- read.csv("../data/classic_fpkm/sample.key.txt",header=T)


## rename fpkm table columns based on sample type
names(tab) <- sample.key$sample_name[match(names(tab), sample.key$sample_id)]
names(tab)[1] <- "tracking_id"

## filter out poorly expressed genes
use <- rowSums(tab[,-1] > 5) >= 3
tab.use <- tab[use,]

## do a log10(FPKM+1) transformation of the data
tab.use[,-1] <- log10(tab.use[,-1]+1)
