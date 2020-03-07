library(WGCNA)
library(randomForestSRC)
library(survival)

ExpressionMatrix <- data.frame(read.table('~/rnaseq.csv', header = T, sep = ','))
rownames(ExpressionMatrix) <- ExpressionMatrix[,1]
ExpressionMatrix <- ExpressionMatrix[,-c(1,2)]
colnames(ExpressionMatrix) <- gsub(".", "-",colnames(ExpressionMatrix), fixed =T)
colnames(ExpressionMatrix) <- gsub("X", "",colnames(ExpressionMatrix), fixed =T)

ClinicalData <- data.frame(read.table("~/clinical_numerical.csv", header = T, sep = ','))
rownames(ClinicalData) <- ClinicalData[,1]
ClinicalData <- ClinicalData[,-1]

Survival <- data.frame(read.table('~/response.csv', header = T, sep = ',', row.names = 1))
ExpressionMatrix <- ExpressionMatrix[, colnames(ExpressionMatrix )%in% rownames(Survival)]
Survival <- Survival[ colnames(ExpressionMatrix ),]
Survival$vitalStatus <- as.numeric(factor(Survival$vitalStatus))
for (i in 1:dim(Survival)[1]) {
  if(Survival$vitalStatus[i] == "2"){
    Survival$vitalStatus[i]<- 1
  } else {
    Survival$vitalStatus[i]<- 0
  }
}
Survival <- as.data.frame(t(Survival))
colnames(Survival) <- colnames(ExpressionMatrix)

# GeneSet.Metzeleret <- as.data.frame(read.csv("~/Desktop/DREAMChallenge_BeatAML/Survival/SurvivalAssociatedGenes_MetzeleretAl.csv", sep=";"))
# GeneSet.Metzeleret <- as.data.frame(GeneSet.Metzeleret[,1])
# GeneSet.Nagyet1 <-as.data.frame(read.csv("~/Desktop/DREAMChallenge_BeatAML/Survival/SurvivalAssociatedgenes_NagyetAl.csv", sep=";"))
# GeneSet.Nagyet1 <- as.data.frame(GeneSet.Nagyet1[!(GeneSet.Nagyet1[,1]%in% ""),  1])
# GeneSet.Nagyet2 <-as.data.frame(read.csv("~/Desktop/DREAMChallenge_BeatAML/Survival/SurvivalAssociatedgenes_NagyetAl_2.csv", sep=";"))
# GeneSet.Nagyet2 <- as.data.frame(GeneSet.Nagyet2[!(GeneSet.Nagyet2[,1]%in% ""),  1])

ExpressionMatrix <- ExpressionMatrix[rowSums(ExpressionMatrix) >5,]

## WGCNA

datExpr0 = as.data.frame(t(ExpressionMatrix))
gsg = goodSamplesGenes(datExpr0, verbose = 3)

# Exploratory Analysis
if (!gsg$allOK){# Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average")

par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)

datExpr0 <- datExpr0[!(rownames(datExpr0) == "14-00800"), ]
ClinicalData <- ClinicalData[rownames(ClinicalData) %in% rownames(datExpr0), ]
Survival <- Survival[, colnames(Survival) != "14-00800"]
sampleTree = hclust(dist(datExpr0), method = "average")

par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)

traitColors = numbers2colors(ClinicalData, signed =F)

plotDendroAndColors(sampleTree, traitColors,groupLabels = names(ClinicalData),main = "Sample dendrogram and trait heatmap")

# Co-correlatio network construction 

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", 
     type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Network construction

net = blockwiseModules(datExpr0, power = 8,TOMType = "unsigned", minModuleSize = 15,reassignThreshold = 0, mergeCutHeight = 0.25, maxBlockSize = 500,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,
                       saveTOMFileBase = "/home/xiaoxiao/Desktop/BeatAML_Challenge/dream_aml/Survival_Fabio/CoExpressionNetwork_BeatAML",verbose = 3)

table(net$colors)

# Plot resulting modules

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)

# Module assignemnt

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
geneTree = net$dendrograms[[1]];
# Eigenvalues for subsequent survival analysis

MEs = net$MEs
save(MEs, moduleLabels, moduleColors, geneTree,file = "/home/xiaoxiao/Desktop/BeatAML_Challenge/dream_aml/Survival_Fabio/CoExpressionNetwork_BeatAML_networkConstruction.RData")

## Correlation module-phenotypical trait

moduleTraitCor = cor(MEs, cbind(ClinicalData, t(Survival[2,])), use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datExpr0))                 

# Visualization correlation

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(cbind(ClinicalData, t(Survival[2,]))),yLabels = names(MEs),ySymbols = names(MEs),colorLabels = FALSE,
               colors = greenWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,
               cex.text = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))

