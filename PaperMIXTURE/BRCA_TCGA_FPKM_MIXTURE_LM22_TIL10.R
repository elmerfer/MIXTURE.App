#We test here the TCGA with MIXER and CIBERSORT.
#
rm(list=ls())
library(edgeR)
library(gplots)
library(data.table)
# library(ComplexHeatmap)
# library(ade4)#distancia de jaccard
library(ggplot2)
# library(circlize)
# library(dtangle)
library(survival)
library(survminer)
library(parallel)
library(openxlsx)
library(immunedeconv)


 library(TCGAbiolinks)
 library(SummarizedExperiment)
# 
library(stringr)

source('Utils/MIXTURE.DEBUG_V0.1.R')
load("Data/LM22.RData")
load("Data/TIL10.RData")
##Prueba
##change the directory to your own directory!!!
#the BRCA RNAseq data can be downloaded from https://www.dropbox.com/s/zki1gkx5mq1quah/BRCA_rna.rds?dl=0
brca.cpms <- readRDS("/home/elmer/Dropbox/Doctorandos/DarioRocha/BRCA/processed_data/BRCA_rna.rds")
TNBC <- apply(brca.cpms$targets[,c("er","pgr","her2")],1, FUN = function(x) all(x == "negative"))
brca.cpms$targets$TNBC <- TNBC ##defining Triple Negative BRCA


##Download FPKM data with TCGAbioLinks
# query_all_brca_fpkm <- GDCquery(project = "TCGA-BRCA",
#                             data.category = "Transcriptome Profiling",
#                             data.type = "Gene Expression Quantification", 
#                             experimental.strategy = "RNA-Seq",
#                             # sample.type = "Primary solid Tumor",
#                             workflow.type = "HTSeq - FPKM")
# 
#  GDCdownload(query_brca_fpkm)
#  GDCdownload(query_all_brca_fpkm)
#  data_all_brca_fpkm <- GDCprepare(query_all_brca_fpkm, save = TRUE, save.filename = "BRCA.All.FPKM.rda")

##Chage the directory according to wahre you save your downloaded data
#path.2.your.data <- "/home/elmer/Dropbox/IDEAS/cibersort/MIXTURE"
# path.2.your.data <- "path/2/my/file"
load(file.path(path.2.your.data,"BRCA.All.FPKM.rda"))
 
 # dim(data)
 # 56499  1222
 gene_annot_all <- rowData(data)
 exp.all <- assay(data)
 target.data.all <- as.data.frame(colData(data))
 tumor.type <- do.call(rbind,str_split(as.character(target.data.all$barcode),"-"))
 table(tumor.type[,4])
 # 01A  01B  01C  06A  11A  11B 
 # 1077   24    1    7   99   14 
 
 dim(gene_annot_all)
 dim(exp.all)
 rownames(exp.all)[1:10]
 tail(gene_annot_all$ensembl_gene_id)
 exp.all <- exp.all[gene_annot_all$ensembl_gene_id,]
 dim(exp.all)
 rownames(exp.all) <- gene_annot_all$external_gene_name
 a.data <- new("EList",list(E = exp.all, genes = gene_annot_all, targets = target.data.all))
 # saveRDS(a.data, file = "Data/BRCA.FPKM.ALL.RDS")
 a.data$targets$TissueType <- do.call(rbind,str_split(as.character(target.data.all$barcode),"-"))[,4]
 table(a.data$targets$TissueType)
 # 01A  01B  01C  06A  11A  11B 
 # 1077   24    1    7   99   14 
 
 
 ##load PBCMC results
 PAM50.PBCMC.results <- readRDS("Data/PAM50.PBCMC.RDS")
 ##Explore and Add PAM50 through PBCMC
dim(PAM50.PBCMC.results)
# [1] 1222    3
pam50.cat <-  do.call(rbind,str_split(PAM50.PBCMC.results$participant,"-"))
table(pam50.cat[,4])
# 01A  01B  01C  06A  11A  11B 
# 1077   24    1    7   99   14

table(PAM50.PBCMC.results$pam50.fpkm.pbcmc)
# Basal   Her2   LumA   LumB    NoA Normal 
# 222    106    236    161    321    176 
 round(100*table(PAM50.PBCMC.results$pam50.fpkm.pbcmc)/sum(table(PAM50.PBCMC.results$pam50.fpkm.pbcmc)),3)
 # Basal   Her2   LumA   LumB    NoA Normal 
 # 18.167  8.674 19.313 13.175 26.268 14.403 

 length(intersect(as.character(target.data.all$barcode), PAM50.PBCMC.results$participant))
 # 1222 !!! OK
 
 a.data$targets <- merge(a.data$targets, PAM50.PBCMC.results[,1:2], by.x = "barcode", by.y = "participant", all.x=T, sort=FALSE)
 
 sum((colnames(a.data$E) == a.data$targets$barcode))
 #1222
 
 # dim(brca.cpms$targets)
 # brca.cpms$targets$barcode <- str_replace_all(brca.cpms$targets$barcode,"\\.","-")
 ##identifico los 1215 que tenemos en la base de cpms
 #hay 1205 en comun
 
 
brca.fpkm <- a.data
table(brca.fpkm$targets$pam50.fpkm.pbcmc)


M.brca.n <- brca.fpkm$E
# df.brca.to.save <- data.frame(genes = rownames(M.brca.n)[rownames(M.brca.n) %in% rownames(LM22)], M.brca.n[rownames(M.brca.n) %in% rownames(LM22),])
# colnames(df.brca.to.save)[1] <- "Gene symbol"
# write.table(df.brca.to.save, file = "/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/BRCA.TCGA.FPKM.txt", quote = F, row.names = F, sep= "\t")
#write.xlsx(df.brca.to.save, file = "/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/BRCA.TCGA.xlsx")

ncores2use <- 10L
######  TO RUN TEST -- Uncomment this section ####
# 
 # brca.cib.n <- MIXTURE(expressionMatrix = M.brca.n, signatureMatrix =  LM22, functionMixture =  cibersort, useCores = ncores2use, verbose  =  TRUE,
 #                       iter = 1000, nullDist = "PopulationBased")
# 
# brca.mix.n <- MIXTURE(expressionMatrix = M.brca.n, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = ncores2use, verbose  =  TRUE,
#                        iter = 1, nullDist = "none")
# TIL10.brca.mix.n <- MIXTURE(expressionMatrix = M.brca.n, signatureMatrix =  TIL10, functionMixture =  nu.svm.robust.RFE, useCores = ncores2use, verbose  =  TRUE,
#                       iter = 1000, nullDist = "none")
# brca.quanti <-deconvolute_quantiseq.default(mix.mat = M.brca.n, 
#                                                  arrays = FALSE, 
#                                                  signame = "/home/elmer/Dropbox/IDEAS/cibersort/MIXTURE/Data/LM22",  
#                                                  tumor = FALSE, 
#                                                  mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
# 
#  TIL10.brca.quanti <-deconvolute_quantiseq.default(mix.mat = M.brca.n, 
#                                                   arrays = FALSE, 
#                                                   signame = "TIL10", 
#                                                   tumor = FALSE, 
#                                                   mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
# TIL10.brca.quanti2 <-deconvolute_quantiseq.default(mix.mat = M.brca.n, 
#                                                   arrays = FALSE, 
#                                                   signame = "TIL10", 
#                                                   tumor = FALSE, 
#                                                   mRNAscale = TRUE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")

# 
# 

# save(brca.mix.n, TIL10.brca.mix.n, brca.quanti,TIL10.brca.quanti,TIL10.brca.quanti2, file ="Data/BRCA.res.RData" ) 
# save(brca.mix.n, TIL10.brca.mix.n, brca.quanti,TIL10.brca.quanti,TIL10.brca.quanti2, file ="Data/BRCA_FPKM.res.1222.RData" ) 
######  (TO RUN TEST) -- Uncomment Up to here ####
load( file ="Data/BRCA_FPKM.res.1222.RData" ) 
#CPMs 
# brca.ciber.web <- ReadCibersortWebResults("Data/BRCA_CIBERSORT.csv", type = "csv")
# TIL10.brca.ciber.web <- ReadCibersortWebResults("Data/CIBERSORT.Output_BRCA.TIL10.csv", type = "csv", nct = 10)
#FPKMs
brca.ciber.web <- ReadCibersortWebResults("Data/CIBERSORT.BRCA.FPKM.1222.LM22.csv", type = "csv")
TIL10.brca.ciber.web <- ReadCibersortWebResults("Data/CIBERSORT.BRCA.FPKM.1222.TIL10.csv", type = "csv", nct = 10)



#Number of detected cell types 
df.brca <- data.frame( Nb = c(apply(GetCellTypes(brca.mix.n),1, sum),
                              apply(brca.ciber.web>0,1,sum),
                              apply(data.matrix(brca.quanti[,-c(1,24)])>0,1, sum)),
                       method = factor(rep(c("MIXTURE","CIBERSORT", "QUANTISEQ"),each = ncol(M.brca.n)), 
                                       levels = c("QUANTISEQ","CIBERSORT", "MIXTURE")))


dp <- ddply(df.brca, .(method), summarise, min= min(Nb, na.rm=T), Q1 = quantile(Nb, 0.25),  Q2 = median(Nb, na.rm=T), Q3 = quantile(Nb, 0.75), max = max(Nb, na.rm=T))
rownames(dp) <- dp[,1]
round(dp[,-1],2)
#CPMs
#           min Q1 Q2 Q3 max
# QUANTISEQ   2 10 11 12  15
# CIBERSORT   8 11 12 13  17
# MIXTURE     2  5  7  8  12
#FPKMs (Total sobre 1205 sujetos)
#           min Q1 Q2 Q3 max
# QUANTISEQ   2  5  8 10  16
# CIBERSORT   8 12 13 14  17
# MIXTURE     2  7  8 10  14
# saveRDS(brca.fpkm, file="BRCA.FPKM.1222.RDS")

##Only primary tumors and not duplicated ones
sum(colnames(brca.fpkm$E) == brca.fpkm$targets$barcode)
# 1222 -> OK!!! they are in the same order
table(brca.fpkm$targets$TissueType)
# 01A  01B  01C  06A  11A  11B 
# 1077   24    1    7   99   14
brca.pt <- brca.fpkm[, which(str_detect(brca.fpkm$targets$TissueType,"01")) ]
dim(brca.pt)
# [1] 56499  1102
table(brca.pt$targets$TissueType)
# 01A  01B  01C 
# 1077   24    1 
## Only primary tumors
# remove duplicated subjects
brca.pt <- brca.pt[, !duplicated(brca.pt$targets$patient)]
dim(brca.pt)
#1091
table(brca.pt$targets$TissueType)
# 01A  01B  01C 
# 1068   22    1 

table(brca.pt$targets$pam50.fpkm.pbcmc, useNA = "always")
# Basal   Her2   LumA   LumB    NoA Normal 
# 205    106    232    160    312     76 
round(100*table(brca.pt$targets$pam50.fpkm.pbcmc, useNA = "always")/sum(table(brca.pt$targets$pam50.fpkm.pbcmc, useNA = "always")),1)
# Basal   Her2   LumA   LumB    NoA Normal   <NA> 
#   18.5   11.0   24.0   11.8   25.6    9.0    0.0 

### Analyzing ONLY PBMC assigned samples ####
##Ploting the amount of estimated cell - types per subject for each method
sum(rownames(GetCellTypes(brca.mix.n)[brca.pt$targets$barcode,]) %in% brca.pt$targets$barcode)
rownames(brca.ciber.web) <- str_replace_all(rownames(brca.ciber.web),"\\.","-" )
sum(rownames(brca.ciber.web[brca.pt$targets$barcode,])%in% brca.pt$targets$barcode)
sum(rownames(brca.quanti[brca.pt$targets$barcode,])%in% brca.pt$targets$barcode)


df.brca <- data.frame( Nb = c(apply(GetCellTypes(brca.mix.n)[brca.pt$targets$barcode,],1, sum),
                              apply(brca.ciber.web[brca.pt$targets$barcode,]>0,1,sum),
                              apply(data.matrix(brca.quanti[,-c(1,24)])[brca.pt$targets$barcode,]>0,1, sum)),
                       method = factor(rep(c("MIXTURE","CIBERSORT", "QUANTISEQ"),each = ncol(brca.pt)), 
                                       levels = c("QUANTISEQ","CIBERSORT", "MIXTURE")))


dp <- ddply(df.brca, .(method), summarise, min= min(Nb, na.rm=T), Q1 = quantile(Nb, 0.25),  Q2 = median(Nb, na.rm=T), Q3 = quantile(Nb, 0.75), max = max(Nb, na.rm=T))
rownames(dp) <- dp[,1]
round(dp[,-1],2)
#             min Q1 Q2 Q3 max
# QUANTISEQ   2  5  8 10  16
# CIBERSORT   8 12 13 14  17
# MIXTURE     2  7  8 10  14
ggplot(df.brca, aes(x=method, y=Nb, fill=method)) +
  geom_boxplot(position=position_dodge(1)) + ggtitle("BRCA")

##Mixture cell types
CT.brca.mix <- cbind(GetMixture(brca.mix.n,"proportion"), brca.fpkm$targets)

CT.cib.web <- cbind(brca.ciber.web,  brca.fpkm$targets)
# CT.brca.abbas <- GetMixture(brca.abbas.n,"proportion")
# CT.brca.dt <- brca.dt.n$estimates
# CT.brca.abis <- GetMixture(brca.abis,"proportion")


##Survival by PAM50 intrinsic sutype (defined by PBCMC - Fresno et al. https://doi.org/10.1093/bioinformatics/btw704)
# clasif <- as.factor(brca.fpkm$targets$pam50.fpkm.pbcmc)
 CT.quanti <- cbind(brca.quanti[,-c(1,24)] , brca.fpkm$targets)



##All samples 
round(100*cbind(QUANTISEQ = apply(brca.quanti[,-c(1,24)],2, function(x) sum(x>0,na.rm=T)),
                                  CIBERSORT=apply(brca.ciber.web,2,function(x) sum(x>0)), MIXTURE=apply(GetMixture(brca.mix.n),2,function(x) sum(x>0)))/nrow(brca.ciber.web),2)

#FPKMs
#                               QUANTISEQ CIBERSORT MIXTURE
# B.cells.naive                    23.32     94.19   55.10
# B.cells.memory                   85.73     13.11    5.73
# Plasma.cells                     99.17     84.81   57.43
# T.cells.CD8                      41.83     96.10   77.93
# T.cells.CD4.naive                 0.17      0.08    0.00
# T.cells.CD4.memory.resting       23.57     97.43   43.49
# T.cells.CD4.memory.activated      1.49     59.34   25.48
# T.cells.follicular.helper        26.39     90.79   58.26
# T.cells.regulatory..Tregs.       22.90     78.59   48.80
# T.cells.gamma.delta              18.26     45.73   18.17
# NK.cells.resting                 19.09     21.49    4.65
# NK.cells.activated                0.58     60.00    8.63
# Monocytes                        21.41     42.57   18.34
# Macrophages.M0                   83.65     88.05   80.50
# Macrophages.M1                   82.82     98.34   95.68
# Macrophages.M2                   57.34     99.83   99.50
# Dendritic.cells.resting          53.44     73.94   44.23
# Dendritic.cells.activated        21.66     28.38   14.44
# Mast.cells.resting               64.65     95.27   67.47
# Mast.cells.activated              6.31      8.05    3.49
# Eosinophils                      18.34      7.72    1.08
# Neutrophils                       7.55     33.36    4.90

##Table 3 of the paper
##Tumor samples
id.tumor <- brca.pt$targets$barcode
id.tumor <- id.tumor[brca.pt$targets$pam50.fpkm.pbcmc %in% c("Basal","Her2","LumA","LumB")]
length(id.tumor)
#703
# id.normal <- which(brca$targets$sample == "11")
dim(GetCellTypes(brca.mix.n)[id.tumor,])
##Only on primary tumors
df.brca <- data.frame( Nb = c(apply(GetCellTypes(brca.mix.n)[id.tumor,],1, sum),
                              apply(brca.ciber.web[id.tumor,]>0,1,sum),
                              apply(data.matrix(brca.quanti[id.tumor,-c(1,24)])>0,1, sum)),
                       method = factor(rep(c("MIXTURE","CIBERSORT", "QUANTISEQ"),each = length(id.tumor)), 
                                       levels = c("QUANTISEQ","CIBERSORT", "MIXTURE")))

dp <- ddply(df.brca, .(method), summarise, min= min(Nb, na.rm=T), Q1 = quantile(Nb, 0.25),  Q2 = median(Nb, na.rm=T), Q3 = quantile(Nb, 0.75), max = max(Nb, na.rm=T))
rownames(dp) <- dp[,1]
round(dp[,-1],2)
#             min Q1 Q2 Q3 max
# QUANTISEQ   2  5  8 10  15
# CIBERSORT   8 12 13 14  17
# MIXTURE     2  7  8 10  13


round(100*cbind(QUANTISEQ = apply(brca.quanti[id.tumor,-c(1,24)],2, function(x) sum(x>0,na.rm=T)),
  CIBERSORT=apply(brca.ciber.web[id.tumor,],2,function(x) sum(x>0)), 
                MIXTURE=apply(GetMixture(brca.mix.n)[id.tumor,],2,function(x) sum(x>0)))/nrow(brca.ciber.web[id.tumor,]),2)

#FPKMs (1091 del primery tumor)
# QUANTISEQ CIBERSORT MIXTURE
# B.cells.naive                    24.56     93.49   56.74
# B.cells.memory                   85.33     14.21    6.14
# Plasma.cells                     98.99     84.05   55.82
# T.cells.CD8                      42.25     95.88   77.54
# T.cells.CD4.naive                 0.18      0.09    0.00
# T.cells.CD4.memory.resting       24.75     97.25   40.60
# T.cells.CD4.memory.activated      1.74     63.06   27.86
# T.cells.follicular.helper        25.57     92.39   61.32
# T.cells.regulatory..Tregs.       25.02     84.05   53.35
# T.cells.gamma.delta              19.89     46.56   16.87
# NK.cells.resting                 18.61     20.71    3.67
# NK.cells.activated                0.73     58.11    6.05
# Monocytes                        19.25     39.78   17.51
# Macrophages.M0                   85.06     92.48   85.70
# Macrophages.M1                   86.16     98.26   95.69
# Macrophages.M2                   58.02     99.82   99.45
# Dendritic.cells.resting          53.44     73.42   43.72
# Dendritic.cells.activated        21.54     26.40   13.84
# Mast.cells.resting               65.54     95.51   66.45
# Mast.cells.activated              5.77      5.87    1.65
# Eosinophils                      15.40      7.42    1.01
# Neutrophils                       8.43     33.73    4.12

round(100*cbind(QUANTISEQ = apply(brca.quanti[id.tumor,-c(1,24)],2, function(x) sum(x>0,na.rm=T)),
                CIBERSORT=apply(brca.ciber.web[id.tumor,],2,function(x) sum(x>0)), 
                MIXTURE=apply(GetMixture(brca.mix.n)[id.tumor,],2,function(x) sum(x>0)))/nrow(brca.ciber.web[id.tumor,]),2)
M.mix <- GetMixture(brca.mix.n)[id.tumor,]
M.cib <- brca.ciber.web[id.tumor,]
M.q <- data.matrix(brca.quanti[id.tumor,-c(1,24)])
M.mix[M.mix==0] <- NA
M.cib[M.cib==0] <- NA
M.q[M.q==0] <- NA
ct.names <- c("Bn","Bm","Pcs","CD8T","CD4T","CD4Tr","CD4Ta","Tfh","Tregs","gdt","Nkr","Nka","MM","M0","M1","M2","Dcr","Dca","Mcr","Mca","E","N")
df.bxp <- data.frame(betahat = c(as.numeric(M.mix), as.numeric(M.cib), as.numeric(M.q)),
                     CT = rep(ct.names, each = nrow(M.mix), times = 3),
                     Method = factor(rep(c("MIXTURE","CIBERSORT","QUANTISEQ"),each=prod(dim(M.mix))),
                                     levels = c("QUANTISEQ","CIBERSORT", "MIXTURE")))

LM22.bx <- ggplot(df.bxp, aes(x=CT, y=betahat, colour = Method)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs( x = " ",y = "Estimated Proportion")+ ggtitle("LM22") +facet_grid(~Method)                     
print(LM22.bx)

ddply(df.bxp, .(Method), summarise, min= min(betahat, na.rm=T), Q1 = quantile(betahat, 0.25, na.rm=T),  Q2 = median(betahat, na.rm=T), 
       Q3 = quantile(betahat, 0.75, na.rm=T), max = max(betahat, na.rm=T))
# Method          min           Q1          Q2          Q3       max
# 1 QUANTISEQ 2.338410e-08 0.0003855035 0.001022389 0.002885611 0.5921421
# 2 CIBERSORT 1.000000e-04 0.0163939388 0.047137960 0.105821095 0.6782144
# 3   MIXTURE 7.033141e-03 0.0348062917 0.077226685 0.162174792 0.9552471



rownames(TIL10.brca.ciber.web) <- str_replace_all(rownames(TIL10.brca.ciber.web), "\\.","-")
TIL10.M.mix <- GetMixture(TIL10.brca.mix.n)[id.tumor,]
TIL10.M.cib <- TIL10.brca.ciber.web[id.tumor,]
TIL10.M.q <- data.matrix(TIL10.brca.quanti[id.tumor,-c(1,12)])
TIL10.M.mix[TIL10.M.mix==0] <-NA
TIL10.M.cib[TIL10.M.cib==0] <- NA
TIL10.M.q[TIL10.M.q==0] <- NA

colnames(TIL10.M.q)
TIL10.ct.names <- c("B","M1","M2","Mo","N","NK","CD4T","CD8T","Tregs","D")
TIL10.df.bxp <- data.frame(betahat = c(as.numeric(TIL10.M.mix), as.numeric(TIL10.M.cib), as.numeric(TIL10.M.q)),
                     CT = rep(TIL10.ct.names, each = nrow(TIL10.M.mix), times = 3),
                     Method = factor(rep(c("MIXTURE","CIBERSORT","QUANTISEQ"),each=prod(dim(TIL10.M.mix))),
                                     levels = c("QUANTISEQ","CIBERSORT", "MIXTURE")))

save(TIL10.df.bxp, file = "DataForbxPlot_BRCA_TIL10.Figura4_1091.RData")
TIL10.bx<-ggplot(TIL10.df.bxp, aes(x=CT, y=betahat, colour = Method)) + geom_boxplot() + 
  labs( x = "Cell types",y = "Estimated Proportion")+ ggtitle("TIL10")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_grid(~Method)                     


# LM22-original	22 = c("Bn","Bm","Pcs","CD8T","CD4T","CD4Tr","CD4Ta","Tfh","Tregs","gdt","Nkr","Nka","MM","M0","M1","M2","Dcr","Dca","Mcr","Mca","E","N")
# LM22-merge4	4	   = c("B", "B", "B",  "CD8T","CD4T","CD4T", "CD4T", "CD4T","CD4T","CD8T","R",  "R",  "R", "R", "R","R",  "R",  "R",  "R",  "R", "R","R")
# then ctClassType = c("B","B","B","CD8T","CD4T","CD4T","CD4T","CD4T","CD4T","CD8T","R","R","R","R","R","R","R","R","R","R","R","R")
# TIL10 example	   = c("B", "B", "B",  "CD8T","CD4T","CD4T", "CD4T", "CD4T","CD4T","CD8T","R",  "R",  "R", "R", "R","R",  "D",  "D",  "R",  "R", "R","R")
LM22.TIL10	 <-     c("B","B", "B",   "CD8T","CD4T","CD4T","CD4T","CD4T","Tregs","CD8T","Nk","Nk","Mo","R","M1","M2","D","D","R","R","R","N")
unique(LM22.TIL10)
colnames(TIL10) 


TIL10.names <-  c("B","M1","M2","Mo","N","NK","CD4T","CD8T","Tregs","D")

M.mix  <- GetMergedCellTypes(GetMixture(brca.mix.n)[id.tumor,], LM22.TIL10)
M.cib <- GetMergedCellTypes(brca.ciber.web[id.tumor,], LM22.TIL10)
M.q <- GetMergedCellTypes(data.matrix(brca.quanti[id.tumor,-c(1,24)]), LM22.TIL10)
M.mix <- M.mix[,colnames(M.mix)!="R"]
M.cib <- M.cib[,colnames(M.cib)!="R"]
M.q <- M.q[,colnames(M.q)!="R"]
M.mix[M.mix==0] <- NA
M.cib[M.cib==0] <- NA
M.q[M.q==0] <- NA

colnames(M.q)
ct.names <- c("Bn","Bm","Pcs","CD8T","CD4T","CD4Tr","CD4Ta","Tfh","Tregs","gdt","Nkr","Nka","MM","M0","M1","M2","Dcr","Dca","Mcr","Mca","E","N")

df.bxp <- data.frame(betahat = c(as.numeric(M.mix), as.numeric(M.cib), as.numeric(M.q)),
                     CT = rep(colnames(M.q), each = nrow(M.mix), times = 3),
                     Method = factor(rep(c("MIXTURE","CIBERSORT","QUANTISEQ"),each=prod(dim(M.mix))),
                                     levels = c("QUANTISEQ","CIBERSORT", "MIXTURE")))
save(df.bxp,TIL10.df.bxp, file = "DataForbxPlot_BRCA_LM22.Figura3a_1091.RData")
LM22.bx <- ggplot(df.bxp, aes(x=CT, y=betahat, colour = Method)) +geom_boxplot() + ylim(0,1) +
  labs( x = "Cell types",y = "Estimated Proportion")+ ggtitle("LM22")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_grid(~Method)                                        
print(LM22.bx)

pdf("/home/elmer/Dropbox/IDEAS/cibersort/GenomeR/Figuras/BRCA_Comp_TIL10_celltypes_Only.pdf", paper = "a4",
    width = 10, height = 10)
grid.arrange(LM22.bx, TIL10.bx, nrow = 2)
dev.off()

##Evaluamos la correlacion entre firmas y metodos
colnames(M.cib)[5] <- "NK"
colnames(M.mix)[5] <- "NK"
colnames(M.q)[5] <- "NK"
colnames(TIL10.M.mix) <- TIL10.names 
colnames(TIL10.M.cib) <- TIL10.names 
colnames(TIL10.M.q)  <- TIL10.names 
summary(cbind(QUANTISEQ=c(M.q),CIBERSORT=c(M.cib),MIXTURE=c(M.mix)))
# QUANTISEQ       CIBERSORT         MIXTURE      
# Min.   :0.000   Min.   :0.0001   Min.   :0.0070  
# 1st Qu.:0.000   1st Qu.:0.0202   1st Qu.:0.0415  
# Median :0.001   Median :0.0652   Median :0.0923  
# Mean   :0.009   Mean   :0.0911   Mean   :0.1202  
# 3rd Qu.:, suggesting     
# Max.   :0.634   Max.   :0.6094   Max.   :0.8250  
# NA's   :3822    NA's   :1381     NA's   :2969
summary(cbind(QUANTISEQ=c(TIL10.M.q),CIBERSORT=c(TIL10.M.cib),MIXTURE=c(TIL10.M.mix)))      
# QUANTISEQ         CIBERSORT        MIXTURE     
# Min.   :0.00000   Min.   :0.000   Min.   :0.007  
# 1st Qu.:0.00000   1st Qu.:0.034   1st Qu.:0.101  
# Median :0.01217   Median :0.108   Median :0.226  
# Mean   :0.03457   Mean   :0.213   Mean   :0.273  
# 3rd Qu.:0.03470   3rd Qu.:0.342   3rd Qu.:0.412  
# Max.   :0.65299   Max.   :0.929   Max.   :0.979  
#                   NA's   :6565    NA's   :7316  
CC.mat <- do.call(rbind, lapply(colnames(M.mix), function(x) {
  c(cor(M.mix[,x],TIL10.M.mix[,x], use = "pairwise.complete.obs"),
    cor(M.cib[,x],TIL10.M.cib[,x], use = "pairwise.complete.obs"),
    cor(M.q[,x],TIL10.M.q[,x], use = "pairwise.complete.obs"))
}))
rownames(CC.mat) <- colnames(M.mix)
colnames(CC.mat) <- c("MIXTURE","CIBERSORT","QUANTISEQ")
summary(CC.mat)
#CPMs
#           MIXTURE          CIBERSORT          QUANTISEQ          
# Min.   :-0.0002549   Min.   :0.0295   Min.   :-0.19288  
# 1st Qu.: 0.1130936   1st Qu.:0.1311   1st Qu.:-0.04651  
# Median : 0.2513331   Median :0.2231   Median :-0.01046  
# Mean   : 0.2608010   Mean   :0.2442   Mean   : 0.04274  
# 3rd Qu.: 0.3917043   3rd Qu.:0.3312   3rd Qu.: 0.17152  
# Max.   : 0.5135764   Max.   :0.5425   Max.   : 0.35687

#FPKM (1091)
#           MIXTURE           CIBERSORT         QUANTISEQ       
# Min.   :-0.29801   Min.   :0.01312   Min.   :-0.03918  
# 1st Qu.: 0.04518   1st Qu.:0.10739   1st Qu.: 0.08305  
# Median : 0.28785   Median :0.26982   Median : 0.25775  
# Mean   : 0.30081   Mean   :0.26907   Mean   : 0.26027  
# 3rd Qu.: 0.48777   3rd Qu.:0.34925   3rd Qu.: 0.40919  
# Max.   : 0.95146   Max.   :0.63609   Max.   : 0.66773  
#                                     NA's   :1        
## 703 FPKM PBCMC Primary Tumor
# MIXTURE          CIBERSORT          QUANTISEQ       
# Min.   :-0.0435   Min.   :-0.11053   Min.   :-0.06098  
# 1st Qu.: 0.1195   1st Qu.: 0.08443   1st Qu.: 0.05349  
# Median : 0.3844   Median : 0.29582   Median : 0.30513  
# Mean   : 0.4248   Mean   : 0.24837   Mean   : 0.26489  
# 3rd Qu.: 0.6514   3rd Qu.: 0.35996   3rd Qu.: 0.39889  
# Max.   : 1.0000   Max.   : 0.60923   Max.   : 0.68819  
#                                     NA's   :1 


# pdf("/home/elmer/Dropbox/IDEAS/cibersort/GenomeR/Figuras/FigurasFinales/CorrMat_BRCA.pdf", paper = "a4", width = 5, height = 5)
# Heatmap(CC.mat, cluster_rows = FALSE, cluster_columns = FALSE, name = "Correlation")
# dev.off()
# png("/home/elmer/Dropbox/IDEAS/cibersort/GenomeR/Figuras/FigurasFinales/CorrMat_BRCA.png")#, paper = "a4", width = 5, height = 5)
# Heatmap(CC.mat, cluster_rows = FALSE, cluster_columns = FALSE, name = "Correlation")
# dev.off()

round(CC.mat,3)

#703
#       MIXTURE CIBERSORT QUANTISEQ
# B       0.695     0.609     0.305
# CD8T    0.445     0.478     0.688
# CD4T    0.324     0.257     0.024
# Tregs   0.312     0.335     0.053
# NK      0.951     0.109    -0.061
# Mo     -0.011     0.353        NA
# M1     -0.043     0.076     0.515
# M2      0.519     0.362     0.399
# D       0.055     0.015     0.306
# N       1.000    -0.111     0.154

wilcox.test(CC.mat[,"MIXTURE"],CC.mat[,"QUANTISEQ"], paired = T, alternative = "greater")
wilcox.test(CC.mat[,"CIBERSORT"],CC.mat[,"QUANTISEQ"], paired = T, alternative = "greater")
wilcox.test(CC.mat[,"MIXTURE"], CC.mat[,"CIBERSORT"], paired = T,  alternative = "greater")



## ALl methods
#we only have 1085 common subjects

brca.cpms$targets$barcode <- str_replace_all(brca.cpms$targets$barcode,"\\.","-")

##let's include ER information
rownames(brca.pt$targets) <- brca.pt$targets$barcode
sum(brca.pt$targets$barcode %in% brca.cpms$targets$barcode)
#1080
sum(brca.cpms$targets$barcode %in% brca.pt$targets$barcode)
brca.pt$targets$ER <- NA
targ.aux <- brca.pt$targets
rownames(targ.aux)[1:2]
dim(targ.aux)
##1091
targ.aux[1:3,1:3]

targ.aux <- merge(targ.aux, brca.cpms$targets[,c("barcode","er")], by.x="barcode", by.y="barcode",all.x=TRUE,all.y=FALSE, sort=FALSE)
rownames(targ.aux) <- targ.aux$barcode

table(targ.aux$er, useNA = "always")
targ.aux$ER[targ.aux$er == "positive"]  <- "ER+"
targ.aux$ER[targ.aux$er == "negative"]  <- "ER-"
table(targ.aux$ER,targ.aux$er, useNA = "always")
# indeterminate negative positive <NA>
#   ER-              0      233        0    0
# ER+              0        0      796    0
# <NA>             2        0        0   60
sum(rownames(brca.pt$targets) == rownames(targ.aux[brca.pt$targets$barcode,]))
brca.pt$targets <- targ.aux[brca.pt$targets$barcode,]



# cbind(rownames(targets.cpm)[1:10],rownames(targets)[1:10])#OK!!!

