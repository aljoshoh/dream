############################################################################################################
# PIPELINE
############################################################################################################


library(foreign)
# Survival Analysis
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(gridExtra)
library(Hmisc)
library(ggplot2)
library(ggrepel)


#Import clinical response
sav <- read.csv("dream_data/response.csv") # Import response
clin_feat_num <- read.csv("dream_data/clinical_categorical.csv")
clin_feat_filter <- clin_feat_num[clin_feat_num$lab_id %in% sav$lab_id, ] # clinical features
row.names(clin_feat_filter) = clin_feat_filter$lab_id
clin_feat_filter <- clin_feat_filter[as.character(sav$lab_id),]
sav$vitalStatus <- as.numeric(factor(sav$vitalStatus))
row.names(sav) = sav$lab_id; sav <- sav[,-1]
sav_single <- cbind(sav, clin_feat_filter)

# Single variant plot
covariate <- "priorMalignancyType"

km_cms_fit <- survfit(Surv(overallSurvival, vitalStatus) ~ sav_single[,covariate], data=sav_single)
cox <- coxph(Surv(overallSurvival, vitalStatus) ~ sav_single[,covariate], data = sav_single)
summary(cox)
autoplot(km_cms_fit)


sav$strat <- "1" #for no sample stratificastion, otherwise import a covariate
drug <- "allacross"
stratification <- "strat"
survivals <- "overallSurvival"
identifier <- "dream_clin_factor_filter"


# SYSTEMATIC
test <- calc_asso_cox(response = sav, features = clin_feat_filter[,-1], survival = survivals, #bem_cox_or
              strat = stratification, treatment = paste("FOLFIRI+",capitalize(drug),sep=""),
              save = paste("metadata/calc_asso_cox_",drug,"_",stratification,"_",survivals,"_",identifier,".RData",sep=""),
              N = 1, plot_volcano = F)



# Importing object &asso_cox

asso_cox_unstrat <- loadRData(paste("metadata/calc_asso_cox_",drug,"_",stratification,"_",survivals,"_",identifier,".RData",sep=""))

asso_cox <- asso_cox_unstrat[[1]] # if ploting the unstratified version
asso_cox$fdr <- p.adjust(asso_cox$P.Value, method = "BH")
sign <- as.numeric(asso_cox$P.Value) < 0.005 & as.numeric(asso_cox$fdr) < 0.5 & abs(as.numeric(asso_cox$effectsize)) > 0.25
sign[is.na(sign)] <- F
print(paste("Found: ",as.character(length(which(sign))), sep=""))

p1 <- plot_volcano(asso_cox, 
                   Cond = sign, 
                   sign = sign, 
                   annotations = T,
                   REFERENCE = asso_cox_unstrat[[1]],
                   PQ = "fdr"
)
plot_volcano_2(asso_cox, Cond = !sign, annotation.object = NULL, pval = 0.05, eff = 0.25,
               max_p = max(c(na.omit(-log10(asso_cox$fdr[sign])),4.4))+0.2,   #!!!! fdr as cutoff
               max_eff = min(max(c(na.omit(abs(asso_cox$effectsize[sign])),1))+0.1,7),
               name = paste(drug," treatment",sep=""),
               PQ = "fdr",
               plot_fdr_cutoff= T)
ggsave(filename = paste("plots/asso_cox_",drug,"_",stratification,"_",survivals,"_",identifier,".pdf",sep=""), width = 10, height = 10)





### Exclusivities
cox_exclusivity <- exclusivity_bem(bem = bem_cox_or, logic = "exclusivity", min_altered = 1)
bem_cox_or <- combinations_bem(bem = bem_cox, logic = "or", min_altered = 10, filter = cox_exclusivity$names[cox_exclusivity$pval < 0.05 & cox_exclusivity$fdr < 0.6])






