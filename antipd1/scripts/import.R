library(readr)
library(tidyverse)
library(Hmisc)

# import
data1c <- read_delim("data/Chen2016_CTLA4_Melanoma_Nanostring/ICB.Chen2016_Ipilimumab_Melanoma.clinical", delim = "\t")
data1m <- read_delim("data/Chen2016_CTLA4_Melanoma_Nanostring/ICB.Chen2016_Ipilimumab_Melanoma.self_subtract", delim = "\t")
data2c <- read_delim("data/Chen2016_PD1_Melanoma_Nanostring_Ipi.Prog/ICB.Chen2016_PD1_Melanoma_Ipi.Prog.clinical", delim = "\t")
data2m <- read_delim("data/Chen2016_PD1_Melanoma_Nanostring_Ipi.Prog/ICB.Chen2016_PD1_Melanoma_Ipi.Prog.self_subtract", delim = "\t")
data3c <- read_delim("data/Gide2019_PD1_Melanoma_RNASeq/ICB.Gide2019_Pembrolizumab-Nivolumab_Melanoma.clinical", delim = "\t")
data3m <- read_delim("data/Gide2019_PD1_Melanoma_RNASeq/ICB.Gide2019_Pembrolizumab-Nivolumab_Melanoma.self_subtract", delim = "\t")
data4c <- read_delim("data/Gide2019_PD1+CTLA4_Melanoma_RNASeq/ICB.Gide2019_Pembrolizumab-Nivolumab+Ipilimumab_Melanoma.clinical", delim = "\t")
data4m <- read_delim("data/Gide2019_PD1+CTLA4_Melanoma_RNASeq/ICB.Gide2019_Pembrolizumab-Nivolumab+Ipilimumab_Melanoma.self_subtract", delim = "\t")
data5c <- read_delim("data/Hugo2016_PD1_Melanoma_RNASeq/ICB.Hugo2016_Pembrolizumab_Melanoma.clinical", delim = "\t")
data5m <- read_delim("data/Hugo2016_PD1_Melanoma_RNASeq/ICB.Hugo2016_Pembrolizumab_Melanoma.self_subtract", delim = "\t")
data6c <- read_delim("data/Kim2018_PD1_Gastric_RNASeq/ICB.Kim2018_Pembrolizumab_Gastric.clinical", delim = "\t")
data6m <- read_delim("data/Kim2018_PD1_Gastric_RNASeq/ICB.Kim2018_Pembrolizumab_Gastric.self_subtract", delim = "\t")
data7c <- read_delim("data/Chen2016_CTLA4_Melanoma_Nanostring/ICB.Chen2016_Ipilimumab_Melanoma.clinical", delim = "\t")
data7m <- read_delim("data/Chen2016_CTLA4_Melanoma_Nanostring/ICB.Chen2016_Ipilimumab_Melanoma.self_subtract", delim = "\t")
data8c <- read_delim("data/Lauss2017_ACT_Melanoma_RNASeq/ICB.Lauss2017_ACT_Melanoma.clinical", delim = "\t")
data8m <- read_delim("data/Lauss2017_ACT_Melanoma_RNASeq/ICB.Lauss2017_ACT_Melanoma.self_subtract", delim = "\t")
data9c <- read_delim("data/Nathanson2017_CTLA4_Melanoma_RNASeq_Post/ICB.Nathanson2017_Ipilimumab_Melanoma_Post.clinical", delim = "\t")
data9m <- read_delim("data/Nathanson2017_CTLA4_Melanoma_RNASeq_Post/ICB.Nathanson2017_Ipilimumab_Melanoma_Post.self_subtract", delim = "\t")
data10c <- read_delim("data/Nathanson2017_CTLA4_Melanoma_RNASeq_Pre/ICB.Nathanson2017_Ipilimumab_Melanoma_Pre.clinical", delim = "\t")
data10m <- read_delim("data/Nathanson2017_CTLA4_Melanoma_RNASeq_Pre/ICB.Nathanson2017_Ipilimumab_Melanoma_Pre.self_subtract", delim = "\t")
data11c <- read_delim("data/Prat2017_PD1_NSCLC-HNSC-Melanoma_Nanostring/ICB.Prat2017_Pembrolizumab-Nivolumab_NSCLC-HNSC-Melanoma.clinical", delim = "\t")
data11m <- read_delim("data/Prat2017_PD1_NSCLC-HNSC-Melanoma_Nanostring/ICB.Prat2017_Pembrolizumab-Nivolumab_NSCLC-HNSC-Melanoma.self_subtract", delim = "\t")
data12c <- read_delim("data/Riaz2017_PD1_Melanoma_RNASeq_Ipi.Naive/ICB.Riaz2017_Nivolumab_Melanoma_Naive.clinical", delim = "\t")
data12m <- read_delim("data/Riaz2017_PD1_Melanoma_RNASeq_Ipi.Naive/ICB.Riaz2017_Nivolumab_Melanoma_Naive.self_subtract", delim = "\t")
data13c <- read_delim("data/Riaz2017_PD1_Melanoma_RNASeq_Ipi.Prog/ICB.Riaz2017_Nivolumab_Melanoma_Prog.clinical", delim = "\t")
data13m <- read_delim("data/Riaz2017_PD1_Melanoma_RNASeq_Ipi.Prog/ICB.Riaz2017_Nivolumab_Melanoma_Prog.self_subtract", delim = "\t")


datac <- list(data1c, data2c, data3c, data4c, data5c, data6c, data7c, data8c, data9c, data10c, data11c, data12c, data13c)
datam <- list(data1m, data2m, data3m, data4m, data5m, data6m, data7m, data8m, data9m, data10m, data11m, data12m, data13m)
datam_ori <- datam

# clean molecular data
for(i in 1:length(datam)){
  colnames(datam[[i]])[1] <- "Entrez"
  tmp <- datam[[i]]
  tmp <- tmp %>% column_to_rownames("Entrez") %>% t %>% as.data.frame %>% rownames_to_column() %>% as.tibble
  datam[[i]] <- tmp
}
datam <- datam %>% reduce(full_join)
datam_save <- datam

#filter genes
# biomart
# filter genes 
library(biomaRt)
dummydata <- read_csv("data/GRCh37ERCC_refseq105_genes_count.csv")
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_list=colnames(datam)
gene_list_dummy <- dummydata$X1
test=getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "entrezgene_id", values = gene_list, bmHeader = T, mart = mart)
new <- lapply(gene_list, function(x) test$`HGNC symbol`[test$`NCBI gene (formerly Entrezgene) ID` == x][1])
new[[1]] <- "row.name"
colnames(datam) = unlist(new)
datam <- datam[,!is.na(colnames(datam))]
datam <- datam[, colnames(datam) %in% c(gene_list_dummy,"row.name")]
immune_genes <- unlist(lapply(data2m$Entrez, function(x) test$`HGNC symbol`[test$`NCBI gene (formerly Entrezgene) ID` == x]))
#datam <- datam[, colnames(datam) %in% c(immune_genes, "row.name")]

genenas <- apply(datam, 2, function(x) length(which(is.na(x))))
datam <- datam[,genenas <65] #<- no filtering for this here
patientnas <- apply(datam, 1, function(x) length(which(is.na(x))))
datam <- datam[patientnas <10,]
datam <- datam[,!duplicated(colnames(datam))]
datam[] <- datam %>% mutate_all(function(x) Hmisc::impute(x))

# clean clinical data                    
for(i in 1:length(datac)){
  colnames(datac[[i]])[1] <- "Patient"
  datac[[i]]$Patient <- as.character(datac[[i]]$Patient)
  if(any(colnames(datac[[i]]) == "RECIST"))
    datac[[i]]$RECIST <- as.character(datac[[i]]$RECIST)
  tmp <- datac[[i]]
  tmp <- tmp %>% as.tibble
  datac[[i]] <- tmp
}
datac <- datac %>% reduce(full_join)
genenas <- apply(datac, 2, function(x) length(which(is.na(x))))
datac <- datac[,genenas <160]
patientnas <- apply(datac, 1, function(x) length(which(is.na(x))))
datac <- datac[patientnas <10,]


# Harmonize and save
datam <- datam[!duplicated(datam$row.name),]
datam <- datam %>% column_to_rownames("row.name")
datac <- datac[!duplicated(datac$Patient),]
datac <- datac %>% column_to_rownames("Patient")
datam <- datam[intersect(row.names(datac), row.names(datam)),]
datac <- datac[intersect(row.names(datac), row.names(datam)),]


save(datam, file = "metadata/features_symbols_v2.RData")
save(datac, file = "metadata/response_v2.RData")


# explore data for batches
labels <- unlist(lapply(row.names(datam), function(x) which(unlist(lapply(datam_ori, function(y) any(x %in% colnames(y)))))[1]))
datap <- datam
datap$set <- labels

library(umap, lib.loc = "/home/alexander.ohnmacht/R/x86_64-redhat-linux-gnu-library/3.6/")
library(reticulate)
library(ggplot2)
var <- apply(datap, 2, sd)
data_red <- datap
data_red$set <- datap$set
use_python("/storage/groups/cbm01/tools/alexander.ohnmacht/miniconda3/bin/python")
custom.config <- umap.defaults
rna.umapa <- umap::umap(data_red, custom.config, method="umap-learn") #_scaled_5000
data_fig1a <- as.data.frame(rna.umapa$layout)
data_fig1a$set <- factor(datap$set)
umap1 <- ggplot(data_fig1a,aes(x=V1,y=V2, group=set)) +
  geom_point(aes(col = set, shape = set))+
  xlab('UMAP1') +
  ylab('UMAP2')+
  scale_fill_viridis_c(guide = "legend")+
  scale_shape_manual(values = c(((0:21)[-c(19,20,21)]),24,25,23))+
  theme_minimal()+
  theme(legend.position = "right")+
  labs(shape = "cohort", col = "cohort")
###############################################
