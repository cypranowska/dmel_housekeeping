# Title: Visualization of normalized data sets for housekeeping gene identificaiton
# Authors: Caroline Cypranowska, Shivali Baveja
# Date: 2019/09/04

library(dplyr)
library(reshape2)
library(tibble)
library(ggplot2)
library(gtools)

option_list <- list(
  make_option("--project", default = "", type = "character", help = "name for project, e.g. bulk_rnai"),
  make_option("--data_dir", default = "../raw", type = "character", help = "directory containing SummarizedExperiment .Rda file"),
  make_option("--out_dir", default = "", type = "character", help = "path to write output"),
  make_option("--gene_list", default = NULL, type = "character", help = "path to .csv file with flybase gene ids for plotting")
)

opt <- parse_args(OptionParser(option_list = option_list))

project <- opt$project
data_dir <- opt$data_dir
out_dir <- opt$out_dir

###----- 1. Load SummarizedExperiment object
load(file = file.path(data_dir,paste0(project,"_se.Rda")))

# create a vector of pairwise combinations of sample IDs
combos <- as.data.frame(combinations(length(colnames(se)),2,colnames(se)))

# fix the levels for each column
combos$V1 <- factor(combos$V1, levels=union(levels(combos$V1), levels(combos$V2)))
combos$V2 <- factor(combos$V2, levels=union(levels(combos$V1), levels(combos$V2)))

# declare a function that makes the log2 CPM plots
## args: key - a dataframe with pairwise sample combinations, data.use - a data frame
## with normalized expression data for all samples
## output: .pdf file for each scatter plot

plotNormExprCombos <- function(key, data.use, out.dir) {
  for (i in 1:nrow(key)) {
    data.use %>% filter(sample_id == key$V1[i] | sample_id == key$V2[i]) %>%
      dcast(flybase_id ~ sample_id) -> tmp
    # remove the flybase_id column
    tmp <- tmp[,-c(1)]
    ggplot(tmp, aes(x = tmp[,1], y = tmp[,2])) + geom_point() + 
      labs(x = paste0("log2(norm counts +1) ", colnames(tmp)[1]), y = paste0("log2(norm counts+1) ", colnames(tmp)[2])) + theme_classic()
    ggsave(filename = file.path(out.dir,paste0("log2normcounts_",colnames(tmp)[1],"_",colnames(tmp)[2],".pdf")))
  }
}

# create an object with the log normalized counts
tab.use <- assays(se)$logcounts

if (is.null(opt$gene_list)) {
  tab.use <- melt(tab.use)
  plotNormExprCombos(combos, tab.use, out.dir = data_dir)
} else if (file.exists(opt$gene_list)) {
  # create an object with the log normalized counts
  gene.list <- read.csv(opt$gene_list, header = T)
  gene.list <- gene.list$flybase_id
  tab.use <- melt(tab.use[match(gene.list, row.names(tab.use)),])
  colnames(tab.use) <- c("flybase_id", "sample_id", "value")
  plotNormExprCombos(combos, tab.use, out.dir = out_dir)
}