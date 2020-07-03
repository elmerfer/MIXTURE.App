#Figure_1c: Methods performance in recognizing individual cell types using LM22 or TIL10 signature
rm(list=ls())
library(immunedeconv)
library(ComplexHeatmap)
library(circlize)
source('Utils/MIXTURE.DEBUG_V0.1.R')
.debug <- TRUE

#Load LM22 Signature
load(file = "Data/LM22.Rdata")
#Load TIL10 Signature
load(file = "Data/TIL10.Rdata")

#Scenario S1

#local implementation of cibersort 
#out.cib <- MIXTURE(expressionMatrix = LM22, signatureMatrix =  LM22, functionMixture =  cibersort, useCores = 3L)
#The CIBERSORT R version provides different results than it's web site counterpart so we dismissed it

#Obtain the outputs from the five evaluated methods using LM22 signature
LM22.output.abbas <- MIXTURE(expressionMatrix = LM22, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = 3L)
LM22.output.abis <- MIXTURE(expressionMatrix = LM22, signatureMatrix =  LM22, functionMixture =  rlm.abis, useCores = 1L)
LM22.output.quantiseq <-deconvolute_quantiseq.default(mix.mat = LM22, arrays = FALSE, signame = "Data/LM22",tumor = FALSE, mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
#Load CIBERSORT output from web results
LM22.output.cibersort.web <- ReadCibersortWebResults(file="Figure_1/Figure_1c/CIBERSORT.Output_LM22.csv", type = "csv")
LM22.output.mixture <- MIXTURE(expressionMatrix = LM22, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = 3L)

#Obtain the outputs from the five evaluated methods using TIL10 signature
TIL10.output.abbas <- MIXTURE(expressionMatrix = TIL10, signatureMatrix =  TIL10, functionMixture =  ls.rfe.abbas, useCores = 3L)
TIL10.output.abis <- MIXTURE(expressionMatrix = TIL10, signatureMatrix =  TIL10, functionMixture =  rlm.abis, useCores = 1L)
TIL10.output.quantiseq <-deconvolute_quantiseq.default(mix.mat = TIL10, arrays = FALSE, signame = "TIL10", tumor = FALSE, mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
#Load CIBERSORT output from web results
TIL10.output.cibersort.web <-  read.csv("Figure_1/Figure_1c/CIBERSORT.Output_TIL10.csv",row.names = "Input.Sample")
TIL10.output.mixture <- MIXTURE(expressionMatrix = TIL10, signatureMatrix =  TIL10, functionMixture =  nu.svm.robust.RFE, useCores = 3L)

#Amount of estimated cell types using LM22 signature
LM22.abbas.fit <- apply(GetCellTypes(LM22.output.abbas)>0, 1, sum)
LM22.abis.fit <- apply(GetMixture(LM22.output.abis)>0,1,sum)
LM22.quantiseq.fit <- apply(LM22.output.quantiseq[,-c(1,24)]>0,1,sum)
LM22.cibersort.fit <- apply(LM22.output.cibersort.web>0,1, sum)
LM22.mixture.fit <- apply(GetCellTypes(LM22.output.mixture)>0,1, sum)

summary(cbind(ABBAS=LM22.abbas.fit,ABIS=LM22.abis.fit,QUANTISEQ=LM22.quantiseq.fit,CIBERSORT=LM22.cibersort.fit,MIXTURE=LM22.mixture.fit))
#ABBAS            ABIS         QUANTISEQ       CIBERSORT        MIXTURE 
#Min.   :6.000   Min.   :10.00   Min.   :1.000   Min.   :1.000   Min.   :1  
#1st Qu.:6.000   1st Qu.:11.00   1st Qu.:1.000   1st Qu.:2.000   1st Qu.:1  
#Median :7.000   Median :12.00   Median :1.000   Median :2.500   Median :1  
#Mean   :7.091   Mean   :11.91   Mean   :1.182   Mean   :3.227   Mean   :1  
#3rd Qu.:8.000   3rd Qu.:13.00   3rd Qu.:1.000   3rd Qu.:4.750   3rd Qu.:1  
#Max.   :9.000   Max.   :13.00   Max.   :2.000   Max.   :8.000   Max.   :1 

#Amount of estimated cell types using LM22 signature
TIL10.abbas.fit <- apply(GetCellTypes(TIL10.output.abbas)>0, 1, sum)
TIL10.abis.fit <- apply(GetMixture(TIL10.output.abis)>0,1,sum)
TIL10.quantiseq.fit <- apply(TIL10.output.quantiseq[,-c(1,12)]>0,1,sum)
TIL10.cibersort.fit <- apply(TIL10.output.cibersort.web>0,1, sum)
TIL10.mixture.fit <- apply(GetCellTypes(TIL10.output.mixture)>0,1, sum)

summary(cbind(ABBAS=TIL10.abbas.fit,ABIS=TIL10.abis.fit,QUANTISEQ=TIL10.quantiseq.fit,CIBERSORT=TIL10.cibersort.fit,MIXTURE=TIL10.mixture.fit))
#ABBAS           ABIS        QUANTISEQ   CIBERSORT       MIXTURE 
#Min.   :4.00   Min.   :4.00   Min.   :1   Min.   :3.00   Min.   :1  
#1st Qu.:4.25   1st Qu.:4.00   1st Qu.:1   1st Qu.:3.25   1st Qu.:1
#Median :5.00   Median :5.50   Median :1   Median :4.00   Median :1  
#Mean   :5.00   Mean   :5.80   Mean   :1   Mean   :4.70   Mean   :1  
#3rd Qu.:5.75   3rd Qu.:7.75   3rd Qu.:1   3rd Qu.:6.00   3rd Qu.:1  
#Max.   :6.00   Max.   :8.00   Max.   :1   Max.   :8.00   Max.   :1  

##Floating point error in ABBAS
abbas.fpe <- GetMixture(LM22.output.abbas,"prop")
diag(abbas.fpe) <- NA
abbas.fpe[abbas.fpe == 0] <- NA
c(min(as.numeric(abbas.fpe),na.rm=T),max(as.numeric(abbas.fpe),na.rm=T))

##Floating point error in ABBIS
abis.fpe <- GetMixture(LM22.output.abis)
diag(abis.fpe) <- NA
c(min(as.numeric(abis.fpe),na.rm=T),max(as.numeric(abis.fpe),na.rm=T))

##Floating point error in CIBERSORT
LM22.output.cibersort.web.fpe <- LM22.output.cibersort.web
diag(LM22.output.cibersort.web.fpe) <- NA
LM22.output.cibersort.web.fpe[LM22.output.cibersort.web.fpe == 0] <- NA
summary(as.numeric(LM22.output.cibersort.web.fpe))

#Extract estimated proportions for each
#for LM22
LM22_abbas_proportions <- GetMixture(LM22.output.abbas, "proportion")
LM22_abis_proportions <- GetMixture(LM22.output.abis)
LM22_quantiseq_proportions <- data.matrix(LM22.output.quantiseq[,-c(1,24)])
LM22_cibersort_proportions <- LM22.output.cibersort.web
LM22_mixture_proportions <- GetMixture(LM22.output.mixture, "proportion")

#Extract estimated proportions for each
#for TIL10
TIL10_abbas_proportions <- GetMixture(TIL10.output.abbas, "proportion")
TIL10_abis_proportions <- GetMixture(TIL10.output.abis)
TIL10_quantiseq_proportions <- data.matrix(TIL10.output.quantiseq[,-c(1,12)])
TIL10_cibersort_proportions <- data.matrix(TIL10.output.cibersort.web[,1:10])
TIL10_mixture_proportions <- GetMixture(TIL10.output.mixture, "proportion")

#Renaming cell types for graphical 
LM22.cell.types.names <- c("BN","BM","PC","CD8","CD4N","CD4Mr","CD4Ma","FH","Tr","TGD","NKr","NKa","M","M0","M1","M2","Dr","Da","Mr","Ma","E","N")
colnames(LM22_abbas_proportions) <- rownames(LM22_abbas_proportions) <- LM22.cell.types.names
colnames(LM22_abis_proportions) <- rownames(LM22_abis_proportions) <- LM22.cell.types.names
colnames(LM22_quantiseq_proportions) <- rownames(LM22_quantiseq_proportions) <- LM22.cell.types.names
colnames(LM22_cibersort_proportions) <- rownames(LM22_cibersort_proportions) <- LM22.cell.types.names
colnames(LM22_mixture_proportions) <- rownames(LM22_mixture_proportions) <- LM22.cell.types.names

TIL10.cell.types.names <- c("B","CD8","CD4","TGD","NK","Mo","D","Mt","Eo","N")
colnames(TIL10_abbas_proportions) <- rownames(TIL10_abbas_proportions) <- TIL10.cell.types.names
colnames(TIL10_abis_proportions) <- rownames(TIL10_abis_proportions) <- TIL10.cell.types.names
colnames(TIL10_quantiseq_proportions) <- rownames(TIL10_quantiseq_proportions) <- TIL10.cell.types.names
colnames(TIL10_cibersort_proportions) <- rownames(TIL10_cibersort_proportions) <- TIL10.cell.types.names
colnames(TIL10_mixture_proportions) <- rownames(TIL10_mixture_proportions) <- TIL10.cell.types.names

#Prepare matrices for visualization
#Define null coefficients
LM22_abbas_proportions[LM22_abbas_proportions == 0] <- NA
LM22_abis_proportions[LM22_abis_proportions <= 0] <- NA
LM22_quantiseq_proportions[LM22_quantiseq_proportions == 0] <- NA
LM22_cibersort_proportions[LM22_cibersort_proportions == 0] <- NA
LM22_mixture_proportions[LM22_mixture_proportions <= 0] <- NA

##For TIL10
TIL10_abbas_proportions[TIL10_abbas_proportions == 0] <- NA
TIL10_abis_proportions[TIL10_abis_proportions <= 0] <- NA
TIL10_quantiseq_proportions[TIL10_quantiseq_proportions == 0] <- NA
TIL10_cibersort_proportions[TIL10_cibersort_proportions == 0] <- NA
TIL10_mixture_proportions[TIL10_mixture_proportions <= 0] <- NA

#Matrix correction for false detected cell types
#ABIS
LM22_abis_proportions_corrected <- LM22_abis_proportions
LM22_abis_proportions_corrected[LM22_abis_proportions_corrected <= 0] <- NA
diag(LM22_abis_proportions_corrected) <- 0 
diag(LM22_abis_proportions_corrected)  <- max(LM22_abis_proportions_corrected, na.rm=T )

#ABBAS
min(as.numeric(LM22_abbas_proportions), na.rm=TRUE)
LM22_abbas_proportions_corrected <- LM22_abbas_proportions
diag(LM22_abbas_proportions_corrected) <- 0.04

summary(cbind(ABBAS= as.numeric(LM22_abbas_proportions_corrected),ABIS=as.numeric(LM22_abis_proportions_corrected), QUANTI = as.numeric(LM22_quantiseq_proportions),CIBERSORT=as.numeric(LM22_cibersort_proportions), MIXTURE = as.numeric(LM22_mixture_proportions)))

#Heatmaps
LM22_heatmap <- 
  Heatmap(LM22_abbas_proportions_corrected, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "ABBAS",name = "ABBAS",row_title = "LM22",row_title_gp = gpar(fontsize=20), column_title_gp = gpar(fontsize=20),col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = TRUE) +
  Heatmap(LM22_abis_proportions_corrected, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "ABIS",name = "ABIS",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(LM22_quantiseq_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "quanTIseq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(LM22_cibersort_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = FALSE, cluster_columns = FALSE, column_title = "CIBERSORT",name = "CIBERSORT",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE, show_column_names = TRUE) +
  Heatmap(LM22_mixture_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = TRUE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F, show_column_names = TRUE)

#Save plot in PDF 12x6
pdf("Figure_1/Figure_1c/LM22_methods_TICT.pdf",width=12,height=6)
print(LM22_heatmap)
dev.off()

TIL10_heatmap <-
  Heatmap(TIL10_abbas_proportions, column_title = "ABBAS", cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE,row_title = "TIL10",row_title_gp = gpar(fontsize=20)) +
  Heatmap(TIL10_abis_proportions, column_title = "ABIS", cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(TIL10_quantiseq_proportions, column_title = "quanTIseq", cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(TIL10_cibersort_proportions, column_title = "CIBERSORT", cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(TIL10_mixture_proportions, column_title = "MIXTURE", cluster_rows = FALSE, show_column_names = TRUE, show_row_names = TRUE, cluster_columns = FALSE,col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = TRUE)

#Save plot in PDF 12x6
pdf("Figure_1/Figure_1c/TIL10_methods_TICT.pdf",width=12,height=6)
print(TIL10_heatmap)
dev.off()



