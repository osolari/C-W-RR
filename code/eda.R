setwd("~/workspace/data/GWAS/code/")
require("ggplot2")
require("dplyr")
require("gplots")
source("./eda.lib.R")
require(corrplot)
require("reshape2")
require("caret")

phenoTypes <- read.table(file = "../variants_omid/spatial_variants/matchedSpatialPhenotypes.csv",
                header = TRUE, sep = "\t")

#pairs(phenoTypes, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel = panel.hist, )


#my_palette <- colorRampPalette(c("blue","black","yellow"))(n = 299)
#heatmap.2(corMat, margin = c(6,6), dendrogram = "row",
 #         col = my_palette, trace="none", density.info = "none", key.title = NA, key.xlab = NA, keysize = .75)

## PCA Plot of Phenotypes

phenoTypesTransform <- preProcess(phenoTypes, method=c("BoxCox", "center", "scale", "pca"))
phenoTypesPCA <- predict(phenoTypesTransform, phenoTypes)
phenoTypesPCA$fish <- "NoFish"
phenoTypesPCA[grepl("-F",rownames(phenoTypesPCA)),"fish"] <- "fish"

outlierSamples <- rownames(phenoTypesPCA[c(order(phenoTypesPCA$PC2, decreasing = T)[c(1:3,length(phenoTypesPCA$PC1))]
                , order(phenoTypesPCA$PC1, decreasing = F)[c(1,length(phenoTypesPCA$PC1))]), ])

phenoTypesPCA$outlier <- NA
phenoTypesPCA[outlierSamples, "outlier"] <- outlierSamples


phenoTypesPCA$site <- unlist(lapply(rownames(phenoTypesPCA), FUN = function(x){return(unlist(strsplit(x, "-"))[1])}))

pdf("./pca.pdf", width = 10, height = 10)
ggplot(data = phenoTypesPCA) + geom_point(aes(x = PC1, y = PC2, col= site, shape = fish)) + 
  geom_text(aes(x = PC1, y = PC2, label = outlier)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank())
dev.off()
  
## Corrplot
phenoTypeFilt <- phenoTypes[!(rownames(phenoTypes) %in% outlierSamples),]
phenoTypesFiltTransObj <- preProcess(phenoTypeFilt, method=c("BoxCox", "center", "scale"))
phenoTypesFiltTrans <- predict(phenoTypesFiltTransObj, phenoTypeFilt)


corMat <- cor(phenoTypesFiltTrans)
pdf("corrplotPhenotype.pdf", width = 10, height = 10)
corrplot(corMat, method = "number", order = "hclust")
#corrplot.mixed(corMat )
dev.off()


## Boxplot of all phenotypes for fish/no-fish condition

phenoTypesFiltTrans$fish <- "NoFish"
phenoTypesFiltTrans[grepl("-F",rownames(phenoTypesFiltTrans)),"fish"] <- "fish"

phenoTypesFiltTrans.m <- melt(phenoTypesFiltTrans, id.vars = "fish")


pdf("boxplot.pdf", width = 15, height = 10)
ggplot(data = phenoTypesFiltTrans.m) + geom_boxplot(aes(x = variable, y = value, fill = fish)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text = element_text(angle = 90, size = 20)) +
  scale_fill_manual(values = c("red","grey")) +
  labs(x = "Phenotypes", y = "Value")
  
dev.off()


## Two-sided Hypothesis testing of the difference in means of different phenotypes

tTest <- phenoTypeTest(data = phenoTypesFiltTrans)
sigPhenotype <- rownames(tTest)[tTest$pVal<0.06]

phenoTypesFiltTransSig <- phenoTypesFiltTrans[,sigPhenotype]

write.table(x = phenoTypesFiltTransSig, file = "~/workspace/data/GWAS/variants_omid/spatial_variants/matchedSpatialPhenotypesFiltTransSig.csv",
            row.names = TRUE, col.names = TRUE, sep = "\t")

X <- read.table(file = "../variants_omid/spatial_variants/matchedSpatialVariants.csv",
                header = TRUE, sep = "\t")
X[!(X == 0)] <- 1
XtransObj <- preProcess(X, method=c("BoxCox", "center", "scale"))
X <- predict(XtransObj, X)
write.table(x = X, file = "~/workspace/data/GWAS/variants_omid/spatial_variants/matchedSpatialVariantsCenteredScaled.csv",
            row.names = TRUE, col.names = TRUE, sep = "\t")

