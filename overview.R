library(org.Hs.eg.db)
library(dplyr)
library(ComplexHeatmap)

library(synapser)
library(synapserutils)
synLogin("aljoshoh","pw")

hs <- org.Hs.eg.db
ensemble <- 
  AnnotationDbi::select(hs, 
       keys = as.character(target_list$target),
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL") # TEC gene is doubled

rnaseq <- read.csv("dream_data/rnaseq.csv")
rna <- rnaseq
row.names(rna) = rnaseq$Gene
rna <- rna[,-c(1,2)] %>% t
dnaseq <- read.csv("dream_data/dnaseq.csv")
auc <- read.csv("dream_data/aucs.csv")
auc <- spread(data=auc, key=inhibitor, value=auc)
row.names(auc) = make.names(auc$lab_id)
auc <- auc[,-1]

clin_response <- read.csv("dream_data/response.csv")
clin_feat_cat <- read.csv("dream_data/clinical_categorical.csv")
clin_feat_legend <- read.csv("dream_data/clinical_categorical_legend.csv")
clin_feat_num <- read.csv("dream_data/clinical_numerical.csv")

### Use string database to filter for informative genes by network propagation of drug targets
### Use morgan fingerprints of drugs that are used



View(list$score)

