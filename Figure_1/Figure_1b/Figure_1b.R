'''Set your working directory in the main folder of GitHub Zip Download'''
#setwd("~/Downloads/MIXTURE.App-master")

#Figure_1b: Auto-correlation matrices
rm(list=ls())
library(ComplexHeatmap)
library(gridExtra)
library(circlize)

#Load LM22 Signature
load(file="Data/LM22.Rdata")
#Load TIL10 Signature
load(file="Data/TIL10.Rdata")

#Calculate autocorrelation of cell-types from the LM22 signature
Correlation_LM22 <- cor(LM22)
#Renaming cell types for graphical purposes
names_LM22 <- c("BN","BM","PC","CD8","CD4N","CD4Mr","CD4Ma","FH","Tr","TGD","NKr","NKa","M","M0","M1","M2","Dr","Da","Mr","Ma","E","N")
colnames(Correlation_LM22) <- names_LM22
rownames(Correlation_LM22) <- names_LM22

#Calculate autocorrelation of cell-types from the TIL10 signature
Correlation_TIL10 <- cor(TIL10)
#Renaming cell types for graphical purposes
names_TIL10 <- c("B","M1","M2","M","N","NK","CD4","CD8","Tr","D")
colnames(Correlation_TIL10) <- names_TIL10
rownames(Correlation_TIL10) <- names_TIL10

#Heatmaps
AC_LM22 <- Heatmap(Correlation_LM22, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = TRUE, cluster_columns = FALSE, 
                    name = paste("Correlation","\nCoefficient"), col = colorRamp2(colors=c("blue", "red"),breaks=c(0,1)),
                    row_names_gp = gpar(fontsize=15,fontfamily="sans"),column_names_gp =gpar(fontsize=15,fontfamily="sans"),
                    column_title_gp = gpar(fontsize=32,fontfamily="sans",fontface=1),column_title = "LM22")
AC_TIL10 <- Heatmap(Correlation_TIL10, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = TRUE, cluster_columns = FALSE, 
                    name = paste("Correlation","\nCoefficient"), col = colorRamp2(colors=c("blue", "red"),breaks=c(0,1)),
                    row_names_gp = gpar(fontsize=15,fontfamily="sans"),column_names_gp =gpar(fontsize=15,fontfamily="sans"),
                    column_title_gp = gpar(fontsize=32,fontfamily="sans",fontface=1),column_title = "TIL10")

#Save PNG 900x600
png(filename = "Figure_1/Figure_1b/Autocorrelation_LM22.png",width = 900,height = 600,units = "px")
AC_LM22
dev.off()

png(filename = "Figure_1/Figure_1b/Autocorrelation_TIL10.png",width = 900,height = 600,units = "px")
AC_TIL10
dev.off()
