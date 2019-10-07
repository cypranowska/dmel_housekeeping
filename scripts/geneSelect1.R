# Title: Identifying housekeeping gene candidates
# Authors: Caroline Cypranowska
# Date: 2019/09/04

library(optparse)
library(dplyr)
library(tibble)
library(SummarizedExperiment)
library(biomaRt)
library(GO.db)
library(ggplot2)
library(GGally)

option_list <- list(
  make_option("--project", default = "", type = "character", help = "name of project, e.g. hs"),
  make_option("--data_dir", default = "../raw", type = "character", help = "path to directory containing SummarizedExperiment .Rda")
)

opt <- parse_args(OptionParser(option_list = option_list))

project <- opt$project
data_dir <- opt$data_dir
out_dir <- "../viz"

###----- 1. Load SummarizedExperiment 
load(file = file.path(data_dir,paste0(project,"_se.Rda")))

###----- 2. Identify candidate genes

# remove genes with any zero raw counts from data set
data.use <- assays(se)$counts
message(paste0("Initial number of genes in data set: ", nrow(data.use)))

keep <- apply(data.use, 1, function(r) all(r != 0))
data.use <- assays(se)$logcounts[keep,]
message(paste0("Genes remaining after removing dimensions with zero values: ", nrow(data.use)))

df <- data.frame(flybase_id = row.names(data.use), 
                 mean.expr = rowMeans(data.use),
                 var.expr = rowVars(data.use))

# filtering genes for against lower bound thresholds for expression
# and upper bound variance

df %>% filter(var.expr < 2.2) %>% filter(mean.expr > 10) -> hk

hk %>% arrange(desc(mean.expr)) %>% arrange(var.expr) %>% top_n(600, -var.expr) -> hk

# plot pairwise comparisons of normalized gene expression
ggpairs(as.data.frame(data.use[which(rownames(data.use) %in% hk$flybase_id),]))
ggsave(file.path(out_dir,"housekeeping_ggpairs_normplot.svg"))

df2 <- as.data.frame(t(data.use[which(rownames(data.use) %in% hk$flybase_id),]))
df2$sample <- rownames(df2)
df2 <- melt(df2)

df2 %>% mutate(tissue = case_when(sample %in% c("SRR070392", "SRR070393", "SRR111884", "SRR111885", "SRR350962", "SRR350963") ~ "im_disc", 
                                  sample %in% c("SRR070405", "SRR070406") ~ "fat_body", 
                                  sample %in% c("SRR070407", "SRR070425") ~ "saliv", 
                                  sample %in% c("SRR070408", "SRR100268") ~ "dig", 
                                  sample %in% c("SRR070409", "SRR070410") ~ "cns", 
                                  sample %in% c("SRR070426", "SRR100269") ~ "carcass")) -> df2

df2 %>% group_by(tissue, variable) %>% summarize(mean_size = mean(value)) %>% dcast(fomula = tissue ~ variable, value.var = "mean_size") -> tissue
tissue %>% column_to_rownames("tissue") %>% t() %>% as.data.frame() -> tissue

ggpairs(tissue)
ggsave(file.path(out_dir,"housekeeping_ggpairs_tissue_normplot.svg"))

###----- 3. Save housekeeping gene lists
write.csv(hk, file = paste0("../filt/housekeeping_final.csv"))

###----- 4. Qualitative observations of selected genes

# get gene ontology terms for genes in hk
fb <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")

hk.go <- getBM(attributes = c('flybase_gene_id', 'external_gene_name','go_id','name_1006','definition_1006'),
               filters = 'flybase_gene_id',
               values = hk$flybase_id,
               mart = fb)
# 
# hk.goslim <- getBM(attributes = c('flybase_gene_id', 'external_gene_name','goslim_goa_accession','goslim_goa_description'),
#                    filters = 'flybase_gene_id',
#                    values = hk$flybase_id,
#                    mart = fb)

go.df <- data.frame(term = Term(hk.go$go_id), goid = names(Term(hk.go$go_id)))
# goslim.df <- data.frame(term = Term(hk.goslim$goslim_goa_accession), goid = names(Term(hk.goslim$goslim_goa_accession)))

# look up gene ontology category for each GO term
ont <- lapply(go.df$goid, function(x) { if (is.na(x)) { 'na' } else { Ontology(GOTERM[[as.character(x)]]) }})
go.df$ont <- unlist(ont)

# ont <- lapply(goslim.df$goid, function(x) { if (is.na(x)) { 'na' } else { Ontology(GOTERM[[as.character(x)]]) }})
# goslim.df$ont <- unlist(ont)

# match flybase ID and GO terms with ontology categories
hk.go %>% left_join(go.df, by = c("go_id" = 'goid')) %>% unique() -> hk.go
#hk.goslim %>% left_join(goslim.df, by = c("goslim_goa_accession" = 'goid')) %>% unique() -> hk.goslim

# subset data on ontology category, sum up instances of each GO term
hk.go %>% filter(ont == "BP") %>% group_by(go_id) %>% tally() %>% arrange(desc(n)) -> hk.go.bp
as.data.frame(hk.go.bp)[1:12,] -> bp.subset
hk.go %>% filter(ont == "MF") %>% group_by(go_id) %>% tally() %>% arrange(desc(n)) -> hk.go.mf
as.data.frame(hk.go.mf)[1:13,] -> mf.subset
hk.go %>% filter(ont == "CC") %>% group_by(go_id) %>% tally() %>% arrange(desc(n)) -> hk.go.cc
as.data.frame(hk.go.cc)[1:10,] -> cc.subset

# plot bar plot of top GO terms for BP category
ggplot(bp.subset, aes(x = reorder(go_id, -n), y = n)) + 
  geom_bar(stat = "identity", color = 'black') + 
  labs(x = 'biological process', y = 'number of genes') + 
  coord_flip() + theme_classic()

ggsave(filename = file.path(out_dir,paste0(project,"final_GO-terms_BP_barplot.svg")))

# plot bar plot of top GO terms for MF category
ggplot(mf.subset, aes(x = reorder(go_id, -n), y = n)) + 
  geom_bar(stat = "identity", color = 'black') + 
  labs(x = 'molecular function', y = 'number of genes') + 
  coord_flip() + theme_classic()

ggsave(filename = file.path(out_dir,paste0(project,"final_GO-terms_MF_barplot.svg")))

# plot bar plot of top GO terms for CC category
ggplot(cc.subset, aes(x = reorder(go_id, -n), y = n)) + 
  geom_bar(stat = "identity", color = 'black') + 
  labs(x = 'cellular component', y = 'number of genes') + 
  coord_flip() + theme_classic()

ggsave(filename = file.path(out_dir,paste0(project,"final_GO-terms_CC_barplot.svg")))

# # repeat for GOSlim terms
# hk.goslim %>% filter(ont == "BP") %>% group_by(goslim_goa_accession) %>% tally() %>% arrange(desc(n)) -> hk.goslim.bp
# as.data.frame(hk.go.bp)[1:16,] -> bp.subset
# hk.goslim %>% filter(ont == "MF") %>% group_by(goslim_goa_accession) %>% tally() %>% arrange(desc(n)) -> hk.goslim.mf
# as.data.frame(hk.goslim.mf)[1:14,] -> mf.subset
# hk.goslim %>% filter(ont == "CC") %>% group_by(goslim_goa_accession) %>% tally() %>% arrange(desc(n)) -> hk.goslim.cc
# as.data.frame(hk.goslim.cc)[1:14,] -> cc.subset
# 
# # plot bar plot of top GOSlim terms for BP category
# ggplot(bp.subset, aes(x = reorder(go_id, -n), y = n)) + 
#   geom_bar(stat = "identity", color = 'black') + 
#   labs(x = 'biological process', y = 'number of genes') + 
#   coord_flip() + theme_classic()
# 
# ggsave(filename = file.path(out_dir,paste0(project,"_GOSlim-terms_BP_barplot.svg")))
# 
# # plot bar plot of top GOSlim terms for MF category
# ggplot(mf.subset, aes(x = reorder(go_id, -n), y = n)) + 
#   geom_bar(stat = "identity", color = 'black') + 
#   labs(x = 'molecular function', y = 'number of genes') + 
#   coord_flip() + theme_classic()
# 
# ggsave(filename = file.path(out_dir,paste0(project,"_GOSlim-terms_MF_barplot.svg")))
# 
# # plot bar plot of top GOSlim terms for CC category
# ggplot(cc.subset, aes(x = reorder(go_id, -n), y = n)) + 
#   geom_bar(stat = "identity", color = 'black') + 
#   labs(x = 'cellular component', y = 'number of genes') + 
#   coord_flip() + theme_classic()
# 
# ggsave(filename = file.path(out_dir,paste0(project,"_GOSlim-terms_CC_barplot.svg")))
