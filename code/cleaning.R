setwd("~/workspace/data/GWAS/code/")
require("ggplot2")
require("stringdist")
require("dplyr")
source("cleaning.lib.R")
require("gridExtra")
spatialVariantsRAW <- read.table(file = "~/workspace/data/GWAS/variants_omid/spatial_variants/dp_4_70_percent_spatial_variants_10000_duplicate_ind_removed_top_293_scaffold_freebayes.recode.parsedScores.phred90.tdl.csv",
                              header = TRUE, sep = "\t", stringsAsFactors = F, check.names = FALSE, row.names = 1)
head(spatialVariantsRAW)
dim(spatialVariantsRAW)
colnames(spatialVariantsRAW) <- gsub("-", "_", colnames(spatialVariantsRAW))
colnames(spatialVariantsRAW)

spatialVariants <- spatialVariantsRAW %>% rename(Blain_F7 = Blaln_F7, Blain_F15 =  Blaln_F15)

phenoType <- read.table(file = "~/workspace/data/GWAS/variants_omid/spatial_variants/spatialVariants.csv"
                        , header = TRUE, check.names = FALSE, sep = ",", stringsAsFactors = FALSE)

rownames(phenoType) <- phenoType$Clone
phenoType <- phenoType %>% dplyr::select(-Pond, -`genome sequenced`, -Clone)
colnames(phenoType) <- c("ageMatP", "ageMother1st", "numOffSpring1st","ageMother2nd", "numOffSpring2nd", "ageMatF1",
                         "sizeMatF1", "ageMother1stF1", "numOffSpring1stF1", "ageMother2ndF1", 
                         "numOffSpring2ndF1", "growth")

matchedNames <- findSimilar(gtList = colnames(spatialVariants), ptList = rownames(phenoType))

spatialVariants <- t(spatialVariants)
X <- spatialVariants[matchedNames$genoType, ]
X <- X[, sort(colSums(X), decreasing = TRUE, index.return = TRUE)$ix[1:10000]]
#DATA <- cbind(matchedSpatialVariants, phenoType)

write.table(x = X, file = "~/workspace/data/GWAS/variants_omid/spatial_variants/matchedSpatialVariants.csv",
            row.names = TRUE, col.names = TRUE, sep = "\t")

write.table(x = phenoType, file = "~/workspace/data/GWAS/variants_omid/spatial_variants/matchedSpatialPhenotypes.csv",
            row.names = TRUE, col.names = TRUE, sep = "\t")

GQ <- read.table(file = "~/workspace/data/GWAS/variants_omid/spatial_variants/dp_4_70_percent_spatial_variants.GQ.tdl.phred90.csv"
                 , header = FALSE)

DP <- read.table(file = "~/workspace/data/GWAS/variants_omid/spatial_variants/dp_4_70_percent_spatial_variants.DP.tdl.phred90.csv"
                 , header = FALSE)

gqdp <- data.frame(GQ, DP)
colnames(gqdp) <- c("GQ", "DP")


p1 <- ggplot(data = gqdp) + geom_histogram(aes(x = GQ, color = DP)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank()) +
  labs(x = "Genotype Quality (Phred Score)", y = "Count")

p2 <- ggplot(data = gqdp) + geom_histogram(aes(x = DP)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank()) +
  labs(x = "Read Depth", y = "Count") + xlim(c(0,30))

p3 <- ggplot(data = gqdp[sample(dim(gqdp)[1], 10000),]) + geom_point(aes(x = DP, y = GQ), alpha = .2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), panel.background = element_blank()) +
  labs(x = "Read Depth", y = "Genotype Quality") + xlim(c(0,30))


pdf("./hist.pdf", width = 10, height = 10)

grid.arrange(p3,p1,p2, ncol = 2, layout_matrix = cbind(c(1,1),c(2,3)))
dev.off()



ggplot(data = gqdp) + geom_histogram(aes(x = GQ, color = DP))

