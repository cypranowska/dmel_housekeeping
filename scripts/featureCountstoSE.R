# Title: featureCounts output to SummarizedExperiment
# Author: Caroline Cypranowska
# Date: 2019/09/02

library(optparse)
library(dplyr)
library(tibble)
library(ggplot2)
library(SummarizedExperiment)

option_list <- list(
  make_option("--project", default = "", type = "character", help = "name for project, e.g. bulk_rnai"),
  make_option("--data_dir", default = "../data/fc_out", type = "character", help = "path to featureCounts output"))

opt <- parse_args(OptionParser(option_list = option_list))

project <- opt$project
data_dir <- opt$data_dir
out_dir <- "../raw"

###----- 1. Load featureCounts output

# first identify all the counts.txt files in data_dir
in_files <- list.files(path = data_dir, pattern = 'counts.txt$', recursive = T)
n.genes <- as.integer(system2("wc", 
                              args = c("-l", file.path(data_dir,in_files[[1]]), " | awk '{print $1}'"),
                              stdout = T)) - 2

# declare a function for loading the counts.txt file
loadfCOut <- function(data.dir, in.file) {
  # read in counts matrix
  data.use <- read.table(file.path(data.dir,in.file), header = T)
  # read in .summary file as column data
  col.data <- as.data.frame((read.table(file.path(data.dir,paste0(in.file,".summary")), header = T)))
  col.data %>% column_to_rownames("Status") %>% t() -> col.data
  rownames(col.data) <- gsub("/counts.txt$","",in.file)
  # take the first 6 columns of data.use as row data
  row.data <- data.use[,c("Geneid","Chr","Start","End","Strand","Length")]
  rownames(row.data) <- row.data$Geneid 
  # keep the last column of data.use
  data.use <- data.use[,ncol(data.use)]
  data.use <- as.matrix(data.use, nrow=length(data.use), ncol=1)
  colnames(data.use) <- rownames(col.data)
  rownames(data.use) <- rownames(row.data)
  # create a temp SummarizedExperiment
  tmp <- SummarizedExperiment(assays = list(counts = data.use), colData = col.data, rowData = row.data)
  return(tmp)
}

# preallocate a Summarized Experiment object for loading in the data
se <- SummarizedExperiment(assays = list(counts = matrix(0,n.genes,length(in_files))))

# iterate over all the sample files
for (i in 1:length(in_files)) {
  tmp <- loadfCOut(data_dir, in_files[[i]])
  assay(se)[,i] <- assay(tmp)
  if (i == 1) {
    rowData(se) <- rowData(tmp)
    col.data <- colData(tmp)
  } else {
    col.data <- rbind(col.data,colData(tmp))
  }
}

# set the row and column names for the SummarizedExperiment object
rownames(se) <- rowData(se)$Geneid
colnames(se) <- gsub("/counts.txt$","",in_files)

colData(se) <- col.data

#####--- 2. Filter SummarizedExperiment 

# declare a CPM function
cpm <- function(dat) {
  y <- apply(dat, 2, function(x) {
    x/sum(x)*1000000
  })
  return(y)
}

# remove genes expressed at < 1 in fewer than 10 cells
keep <- rowSums(cpm(assay(se)) >= 5) > 10
se <- se[keep,]
assays(se)$normcounts <- cpm(assay(se))
assays(se)$logcounts <- log2(cpm(assay(se)+1))

#####--- 3. Save SummarizedExperiment object
save(se, file = file.path(out_dir,paste0(project,"_se.Rda")))
