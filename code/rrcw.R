setwd("~/workspace/data/GWAS/code/")
require("MASS")
require("ggplot2")
require("randomForest")
require("glmnet")
require("dplyr")
source("rrcw.lib.R")
require("reshape2")
require("caret")

X <- read.table(file = "../variants_omid/spatial_variants/matchedSpatialVariantsCenteredScaled.csv",
                header = TRUE, sep = "\t")

Y <- read.table(file = "../variants_omid/spatial_variants/matchedSpatialPhenotypesFiltTransSig.csv",
                header = TRUE, sep = "\t")
rownames(Y) <- tolower(gsub(pattern = "-", replacement = "_", rownames(Y)))
rownames(X) <- tolower(rownames(X))
X <- X[rownames(Y),]

TSEtot <- data.frame()
for (l in 1:50){
  
  print(l)
  IXtrn <- sample(1:dim(X)[1], 100)
  IXtst <- setdiff(1:dim(X)[1], IXtrn)

  Xtrn <- X[IXtrn,]
  Ytrn <- Y[IXtrn,]
  Xtst <- X[IXtst,]
  Ytst <- Y[IXtst,]

  YHAT <- predictAllPhenotypes(x = Xtrn, y = Ytrn, xnew = Xtst, method = "ridge")
  dropCol <- apply(YHAT$Yhattrn, 2, FUN = function(x){length(unique(x))==1})

  Yhattrn <- YHAT$Yhattrn[,!dropCol]
  Yhattst <- YHAT$Yhattst[,!dropCol]
  Ytrn <- Ytrn[,!dropCol]
  Ytst <- Ytst[,!dropCol]

  pN <- dim(Xtrn)[2]/dim(Xtrn)[1]

  tdtGCV <- rrCWGCV(y = Ytrn, yhat = Yhattrn, r = pN )
  
  tdtFCV <- rrCWFCV(x = Xtrn, y = Ytrn, yhat = Yhattrn, nfold = 5)
  
  YtildtstGCV <- Yhattst %*% tdtGCV
  YtildtstFCV <- Yhattst %*% tdtFCV

  rrcwTSEGCV <- data.frame(colSums((YtildtstGCV - Ytst)^2))
  rrcwTSEFCV <- data.frame(colSums((YtildtstFCV - Ytst)^2))
  rrTSE <- data.frame(colSums((Yhattst - Ytst)^2))
  
  simTSE <- cbind(rrTSE, rrcwTSEGCV, rrcwTSEFCV) 
  simTSE$phen <- row.names(simTSE)
  
  TSEtot <- rbind(TSEtot, simTSE)
  
  }


head(TSEtot)
unique(TSEtot$phen)
colnames(TSEtot) <- c("RR", "CWRR_GCV", "CWRR_FCV", "phenotype")

ageMatF1 <- TSEtot[TSEtot$phen == "ageMatF1",]

ggplot(data = melt(ageMatF1)) + geom_boxplot(aes(x = variable, y = value)) 

TSEmelt <- melt(data = TSEtot, id.vars = "phenotype")

pdf("./boxplotMethods.pdf", width = 10, height = 10)
ggplot(data = TSEmelt) + geom_boxplot(aes(x = variable, y = value, fill = variable)) + facet_grid(.~phenotype) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(angle = 90, size = 10)) +
  scale_fill_manual(values = c("red","grey", "white")) +
  labs(x ="", y = "Total Squared Error")
dev.off()

##
## Fish/No-Fish Samples

Yf <- Y[grepl(pattern = "_f", x = rownames(Y)),]
Ynf <- Y[grepl(pattern = "_nf", x = rownames(Y)),]

Y <- Yf
X <- X[rownames(Yf),]
