#Supplementary Figure 2: Methods benchmarking on pure cell lines data
library(openxlsx)
library(ggplot2)
library(MIXTURE)
library(ggpubr)
source("Utils/MIXTURE.DEBUG_V0.1.R")

rm(list=ls())

#Generate plotting function
gplotcells <- function(df,title=NULL){
  ggplot(df, aes(x=CT,y=Y)) + geom_boxplot(outlier.size = 0.25) +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, size=6,hjust = 1),
          strip.background = element_rect(colour="black"),
          axis.title.x = element_blank()) +
    ylab("Estimated cell type proportion") + 
    facet_wrap(~Method, nrow=1) +
    if(is.null(title)==FALSE) {
      ggtitle(label=title)
    }
}
# Generate function to create dataframes
makedf <- function(mat, name, sig){
  data.frame(Y = c(data.matrix(mat)), CT = rep(colnames(mat),each=nrow(mat)),Method = name, Signature = sig)
}

# [1] "B.cells.naive"                "B.cells.memory"               "Plasma.cells"                
# [4] "T.cells.CD8"                  "T.cells.CD4.naive"            "T.cells.CD4.memory.resting"  
# [7] "T.cells.CD4.memory.activated" "T.cells.follicular.helper"    "T.cells.regulatory..Tregs."  
# [10] "T.cells.gamma.delta"          "NK.cells.resting"             "NK.cells.activated"          
# [13] "Monocytes"                    "Macrophages.M0"               "Macrophages.M1"              
# [16] "Macrophages.M2"               "Dendritic.cells.resting"      "Dendritic.cells.activated"   
# [19] "Mast.cells.resting"           "Mast.cells.activated"         "Eosinophils"                 
# [22] "Neutrophils"                
colnames.map.lm22 <- c("BN","BM","PC","CD8","CD4N","CD4Mr","CD4Ma","Tcfh","Tregs","Tcgd","NKr","NKa",
                       "Mo","M0","M1","M2","DR","DA","MA","MR","E","N")

##Cell lines analyzed by ABBAS
abbas.cell.lm22 <- readRDS(file = "Figure_2/Supplementary_Figure2/ABBAS.cells.LM22.RDS")
abbas.cell.til10 <- readRDS(file = "Figure_2/Supplementary_Figure2/ABBAS.cells.TIL10.RDS")

##Cell lines analyzed by ABIS
abis.cell.lm22 <- readRDS(file = "Figure_2/Supplementary_Figure2/ABIS.cells.LM22.RDS")
abis.cell.til10 <- readRDS(file = "Figure_2/Supplementary_Figure2/ABIS.cells.TIL10.RDS")

##Cell lines analyzed by quanTIseq
quanTIseq.cell.lm22 <- readRDS(file = "Figure_2/Supplementary_Figure2/quanTIseq.cells.LM22.RDS")
quanTIseq.cell.til10 <- readRDS(file = "Figure_2/Supplementary_Figure2/quanTIseq.cells.TIL10.RDS")

##Cell lines analyzed by CIBERSORT
cibersort.cell.lm22<- readRDS(file="Figure_2/Supplementary_Figure2/CIBERSORT.cells.LM22.RDS")
cibersort.cell.til10<- readRDS(file="Figure_2/Supplementary_Figure2/CIBERSORT.cells.TIL10.RDS")

colnames(cibersort.cell.lm22)
rownames(cibersort.cell.lm22) <- cibersort.cell.lm22[,1]
cibersort.cell.lm22 <- cibersort.cell.lm22[,-c(1,24:ncol(cibersort.cell.lm22))]

rownames(cibersort.cell.til10) <- cibersort.cell.til10[,1]
cibersort.cell.til10 <- cibersort.cell.til10[,-c(1,12:ncol(cibersort.cell.til10))]
colnames(cibersort.cell.til10)

##Cell lines analyzed by MIXTURE (RFE-v-svr)
mixture.cell.lm22 <- readRDS(file = "Figure_2/Supplementary_Figure2/MIXTURE.cells.LM22.RDS")
mixture.cell.til10 <- readRDS(file = "Figure_2/Supplementary_Figure2/MIXTURE.cells.TIL10.RDS")

##Distribution of absolute values
##LM22
abbas.cell.lm22 <- GetMixture(abbas.cell.lm22,type = "absolute")
abis.cell.lm22 <- GetMixture(abis.cell.lm22,type = "absolute")
quanTIseq.cell.lm22 <- quanTIseq.cell.lm22$Abs[,-ncol(quanTIseq.cell.lm22$Abs)]
cibersort.cell.lm22  <- data.matrix(cibersort.cell.lm22)
mixture.cell.lm22 <- GetMixture(mixture.cell.lm22,type = "absolute")

#Metrics of number of 0s estimated for each cell line by each method
summary(colSums(abbas.cell.lm22==0)/nrow(mixture.cell.lm22))
summary(colSums(abis.cell.lm22==0)/nrow(mixture.cell.lm22))
summary(colSums(quanTIseq.cell.lm22==0)/nrow(mixture.cell.lm22))
summary(colSums(cibersort.cell.lm22==0)/nrow(mixture.cell.lm22))
mixture.cell.lm22[is.na(mixture.cell.lm22)]<-0
summary(colSums(mixture.cell.lm22==0)/nrow(mixture.cell.lm22))

# Set mapped cell types for LM22
colnames(abbas.cell.lm22) <- colnames.map.lm22
colnames(abis.cell.lm22) <- colnames.map.lm22
colnames(quanTIseq.cell.lm22) <- colnames.map.lm22
colnames(cibersort.cell.lm22) <- colnames.map.lm22
colnames(mixture.cell.lm22) <- colnames.map.lm22

# Generate dataframe
dfct.lm22 <- rbind(
  makedf(data.matrix(abbas.cell.lm22), "ABBAS", "LM22"),
  makedf(data.matrix(abis.cell.lm22), "ABIS", "LM22"),
  makedf(quanTIseq.cell.lm22, "quanTIseq", "LM22"),
  makedf(cibersort.cell.lm22, "CIBERSORT", "LM22"),
  makedf(data.matrix(mixture.cell.lm22), "MIXTURE", "LM22")
)

# Plot results for LM22 Signature
LM22.panel <- gplotcells(df = dfct.lm22,title = "LM22")
LM22.panel
dev.off()

# [1] "B.cells"         "Macrophages.M1"  "Macrophages.M2"  "Monocytes"       "Neutrophils"     "NK.cells"       
# [7] "T.cells.CD4"     "T.cells.CD8"     "Tregs"           "Dendritic.cells"
colnames.til10.map <-c("B","M1","M2","Mo","N","NK","CD4","CD8","Tregs","D")

##Distribution of absolute values
##TIL10
abbas.cell.til10 <- GetMixture(abbas.cell.til10,type = "absolute")
abis.cell.til10 <- GetMixture(abis.cell.til10,type = "absolute")
quanTIseq.cell.til10 <- quanTIseq.cell.til10$Abs[,-ncol(quanTIseq.cell.til10$Abs)]
cibersort.cell.til10 <- data.matrix(cibersort.cell.til10)
mixture.cell.til10 <-  GetMixture(mixture.cell.til10,type = "absolute")

#Metrics of number of 0s estimated for each cell line by each method
summary(colSums(abbas.cell.til10==0)/nrow(mixture.cell.til10))
summary(colSums(abis.cell.til10==0)/nrow(mixture.cell.til10))
summary(colSums(quanTIseq.cell.til10==0)/nrow(mixture.cell.til10))
summary(colSums(cibersort.cell.til10==0)/nrow(mixture.cell.til10))
mixture.cell.til10[is.na(mixture.cell.til10)]<-0
summary(colSums(mixture.cell.til10==0)/nrow(mixture.cell.til10))

# Set mapped cell types for TIL10
colnames(abbas.cell.til10) <- colnames.til10.map
colnames(abis.cell.til10) <- colnames.til10.map
colnames(quanTIseq.cell.til10) <- colnames.til10.map
colnames(cibersort.cell.til10) <- colnames.til10.map
colnames(mixture.cell.til10) <- colnames.til10.map

# Generate dataframe
dfct.til10 <- rbind(
  makedf(abbas.cell.til10, "ABBAS", "TIL10"),
  makedf(abis.cell.til10, "ABIS", "TIL10"),
  makedf(quanTIseq.cell.til10, "quanTIseq", "TIL10"),
  makedf(cibersort.cell.til10, "CIBERSORT", "TIL10"),
  makedf(mixture.cell.til10, "MIXTURE", "TIL10")
)

# Round data
dfct.til10$Y <- round(dfct.til10$Y,2)
# Plot results for TIL10 Signature
TIL10.panel <- gplotcells(df = dfct.til10,title = "TIL10")
TIL10.panel
dev.off()

# Save Supplementary Figure 2 in PDF 12X10
ggsave(plot = ggarrange(plotlist = list(LM22.panel,TIL10.panel),ncol = 1,nrow = 2),filename = "Figure_2/Supplementary_Figure2/Supplementary_Figure2.pdf",device = "pdf",width = 12,height = 10)


