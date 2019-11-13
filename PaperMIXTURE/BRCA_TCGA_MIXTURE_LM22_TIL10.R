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
library(dtangle)
library(survival)
library(survminer)
library(parallel)
library(openxlsx)
library(immunedeconv)
library(gridExtra)

source('Utils/MIXTURE.DEBUG_V0.1.R')
load("Data/LM22.RData")
load("Data/TIL10.RData")
##Prueba
##change the directory to your own directory!!!
#the BRCA RNAseq data can be downloaded from https://www.dropbox.com/s/zki1gkx5mq1quah/BRCA_rna.rds?dl=0
brca <- readRDS("/home/elmer/Dropbox/Doctorandos/DarioRocha/BRCA/processed_data/BRCA_rna.rds")
TNBC <- apply(brca$targets[,c("er","pgr","her2")],1, FUN = function(x) all(x == "negative"))
brca$targets$TNBC <- TNBC ##defining Triple Negative BRCA


FindSurvCutoff <- function(df.fs, n.intervals,plot=FALSE,title=""){
  cutoff <- seq(min(df.fs$CT),max(df.fs$CT),length.out = n.intervals)[-c(1,n.intervals)]
  pvals <- unlist(lapply(cutoff, function(i){
    df.fs$Q <- "H"
    df.fs$Q[df.fs$CT <= i] <- "L"
    fit.m <- try(survdiff(Surv(Time, Event) ~ Q, data = df.fs))
    if(class(fit.m) == "try-error") return(NA)
    pchisq(fit.m$chisq, length(fit.m$n)-1, lower.tail = FALSE)
  }))
  
  ic <- cutoff[which.min(pvals)]
  pval <- pvals[which.min(pvals)]
  df.fs$Q <- "H"
  df.fs$Q[df.fs$CT <= ic] <- "L"
  fit.m <- survfit(Surv(Time, Event) ~ Q, data = df.fs)
  gs <- ggsurvplot(fit.m, df.fs, pval = TRUE, risk.table = T) + ggtitle(paste(round(ic,3),title,sep=" - "))
  
  if(plot) print(gs)
  invisible(return(list(gplot=gs, cutoff = ic, pval = pval, surv.summary = surv_summary(fit.m))))
  # p = 0.00226
  # e=0.1
  # z <- qnorm(1-(0.5*p))
  # (4*dnorm(z)/z)+dnorm(z)*(z-(1/z))*log((1-e)*(1-e)/(e*e))
  
}
CoxSurvival <- function(df.fs){

  coxph.summary <- summary(coxph(Surv(Time, Event) ~ CT, data = df.fs))
  return(coxph.summary$coefficients)
  
  # p = 0.00226
  # e=0.1
  # z <- qnorm(1-(0.5*p))
  # (4*dnorm(z)/z)+dnorm(z)*(z-(1/z))*log((1-e)*(1-e)/(e*e))
  
}


quantilize <- function(var, qs = c(0.25,0.5,0.75)){
  varq <- quantile(var,qs,na.rm=T)
  if(all(varq==0)){
    varf <- rep("<med", length(var))
    varf[var > mean(var)] <- ">med"
    return(factor(varf))
  }
  if(varq[1] == 0) varq[1]<-0.00001
  
  varn <- paste("Q", 1:(length(var)+1),sep="")
  varf <- rep(varn[1], length(var))
  for( i in 2:length(varn)){
    varf[var >= varq[i-1]] <- varn[i]
  }
  return(factor(varf))
}

##normalize brca counts
 dge <- DGEList(counts = brca$E)
 
 # dge <- calcNormFactors(dge)
 brca.norm <- brca
 brca.norm$E <- cpm(dge$counts) #library size normalization
# 
 M.brca.n <- brca.norm$E
 sum(rownames(M.brca.n) %in% rownames(LM22))
#df.brca.to.save <- data.frame(genes = rownames(M.brca.n)[rownames(M.brca.n) %in% rownames(LM22)], M.brca.n[rownames(M.brca.n) %in% rownames(LM22),])
#colnames(df.brca.to.save)[1] <- "Gene symbol"
# write.table(df.brca.to.save, file = "/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/BRCA.TCGA.txt", quote = F, row.names = F, sep= "\t")
#write.xlsx(df.brca.to.save, file = "/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/BRCA.TCGA.xlsx")

ncores2use <- 4L
# 
 # brca.cib.n <- MIXTURE(expressionMatrix = M.brca.n, signatureMatrix =  LM22, functionMixture =  cibersort, useCores = ncores2use, verbose  =  TRUE,
 #                       iter = 1000, nullDist = "PopulationBased")
# 
brca.mix.n <- MIXTURE(expressionMatrix = M.brca.n, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = ncores2use, verbose  =  TRUE,
                       iter = 1, nullDist = "none")
TIL10.brca.mix.n <- MIXTURE(expressionMatrix = M.brca.n, signatureMatrix =  TIL10, functionMixture =  nu.svm.robust.RFE, useCores = ncores2use, verbose  =  TRUE,
                      iter = 1000, nullDist = "none")
brca.quanti <-deconvolute_quantiseq.default(mix.mat = M.brca.n, 
                                                 arrays = FALSE, 
                                                 signame = "/home/elmer/Dropbox/IDEAS/cibersort/MIXTURE/Data/LM22",  
                                                 tumor = FALSE, 
                                                 mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")

 TIL10.brca.quanti <-deconvolute_quantiseq.default(mix.mat = M.brca.n, 
                                                  arrays = FALSE, 
                                                  signame = "TIL10", 
                                                  tumor = FALSE, 
                                                  mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
TIL10.brca.quanti2 <-deconvolute_quantiseq.default(mix.mat = M.brca.n, 
                                                  arrays = FALSE, 
                                                  signame = "TIL10", 
                                                  tumor = FALSE, 
                                                  mRNAscale = TRUE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")

# 
# 

# save(brca.mix.n, TIL10.brca.mix.n, brca.quanti,TIL10.brca.quanti,TIL10.brca.quanti2, file ="Data/BRCA.res.RData" ) 
load( file ="Data/BRCA.res.RData" ) 
# 
brca.ciber.web <- ReadCibersortWebResults("Data/BRCA_CIBERSORT.csv", type = "csv")
TIL10.brca.ciber.web <- ReadCibersortWebResults("Data/CIBERSORT.Output_BRCA.TIL10.csv", type = "csv", nct = 10)




#Number of detected cell types 
df.brca <- data.frame( Nb = c(apply(GetCellTypes(brca.mix.n),1, sum),
                              apply(brca.ciber.web>0,1,sum),
                              apply(data.matrix(brca.quanti[,-c(1,24)])>0,1, sum)),
                       method = factor(rep(c("MIXTURE","CIBERSORT", "QUANTISEQ"),each = ncol(M.brca.n)), 
                                       levels = c("QUANTISEQ","CIBERSORT", "MIXTURE")))


dp <- ddply(df.brca, .(method), summarise, min= min(Nb, na.rm=T), Q1 = quantile(Nb, 0.25),  Q2 = median(Nb, na.rm=T), Q3 = quantile(Nb, 0.75), max = max(Nb, na.rm=T))
rownames(dp) <- dp[,1]
round(dp[,-1],2)
#           min Q1 Q2 Q3 max
# QUANTISEQ   2 10 11 12  15
# CIBERSORT   8 11 12 13  17
# MIXTURE     2  5  7  8  12

##Ploting the amount of estimated cell - types per subject for each method
ggplot(df.brca, aes(x=method, y=Nb, fill=method)) +
  geom_boxplot(position=position_dodge(1)) + ggtitle("BRCA")

##Mixture cell types
CT.brca.mix <- cbind(GetMixture(brca.mix.n,"proportion"), brca$targets)
CT.cib.web <- cbind(brca.ciber.web, , brca$targets)
# CT.brca.abbas <- GetMixture(brca.abbas.n,"proportion")
# CT.brca.dt <- brca.dt.n$estimates
# CT.brca.abis <- GetMixture(brca.abis,"proportion")


##Survival by PAM50 intrinsic sutype (defined by PBCMC - Fresno et al. https://doi.org/10.1093/bioinformatics/btw704)
clasif <- as.factor(brca$targets$pam50.pbcmc)
CT.quanti <- cbind(brca.quanti[,-c(1,24)] , brca$targets)



##All samples 
round(100*cbind(QUANTISEQ = apply(brca.quanti[,-c(1,24)],2, function(x) sum(x>0,na.rm=T)),
                                  CIBERSORT=apply(brca.ciber.web,2,function(x) sum(x>0)), MIXTURE=apply(GetMixture(brca.mix.n),2,function(x) sum(x>0)))/nrow(brca.ciber.web),2)
write.table(round(100*cbind(CIBERSORT=apply(brca.ciber.web,2,function(x) sum(x>0)), MIXTURE=apply(GetMixture(brca.mix.n),2,function(x) sum(x>0)))/nrow(brca.ciber.web),2), file ="Tabla.txt")
#                               QUANTISEQ CIBERSORT MIXTURE
# B.cells.naive                    83.62     91.93   35.14
# B.cells.memory                    8.81     17.37    5.84
# Plasma.cells                     93.25     82.88   60.99
# T.cells.CD8                       8.07     72.10   31.93
# T.cells.CD4.naive                 0.16      1.81    0.33
# T.cells.CD4.memory.resting       66.58     96.38   42.39
# T.cells.CD4.memory.activated      1.15     17.28    3.21
# T.cells.follicular.helper        81.98     99.09   90.95
# T.cells.regulatory..Tregs.       78.85     83.70   44.12
# T.cells.gamma.delta               0.91      6.58    0.91
# NK.cells.resting                 32.43     77.37   14.07
# NK.cells.activated                0.33     43.79    3.70
# Monocytes                        57.37     57.86   10.62
# Macrophages.M0                   85.02     83.13   73.66
# Macrophages.M1                   93.66     91.19   74.81
# Macrophages.M2                   92.51     99.26   97.78
# Dendritic.cells.resting          78.85     27.57    7.65
# Dendritic.cells.activated        45.93     36.30   17.86
# Mast.cells.resting               88.48     96.71   50.53
# Mast.cells.activated              4.12      4.20    0.91
# Eosinophils                      50.62      6.67    0.00
# Neutrophils                      14.65     12.18    1.23

#                           CIBERSORT MIXTURE
# B cells naive                91.93 35.14
# B cells memory               17.37  5.84
# Plasma cells                 82.88 60.99 *
# T cells CD8                  72.10 31.93
# T cells CD4 naive             1.81  0.33
# T cells CD4 memory resting   96.38 42.39
# T cells CD4 memory activated 17.28  3.21
# T cells follicular helper    99.09 90.95  *https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3696556/
# T cells regulatory  Tregs    83.70 44.12
# T cells gamma delta           6.58  0.91
# NK cells resting             77.37 14.07
# NK cells activated           43.79  3.70
# Monocytes                    57.86 10.62
# Macrophages M0               83.13 73.66 
# Macrophages M1               91.19 74.81 *https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5509958/
# Macrophages M2               99.26 97.78
# Dendritic cells resting      27.57  7.65
# Dendritic cells activated    36.30 17.86
# Mast cells resting           96.71 50.53 * https://link.springer.com/article/10.1007%2Fs10549-007-9546-3
# Mast cells activated          4.20  0.91
# Eosinophils                   6.67  0.00
# Neutrophils                  12.18  1.23

##Table 3 of the paper
##Tumor samples
id.tumor <- which(brca$targets$sample == "01")
id.normal <- which(brca$targets$sample == "11")
round(100*cbind(QUANTISEQ = apply(brca.quanti[id.tumor,-c(1,24)],2, function(x) sum(x>0,na.rm=T)),
  CIBERSORT=apply(brca.ciber.web[id.tumor,],2,function(x) sum(x>0)), 
                MIXTURE=apply(GetMixture(brca.mix.n)[id.tumor,],2,function(x) sum(x>0)))/nrow(brca.ciber.web[id.tumor,]),2)
write.table(round(100*cbind(CIBERSORT=apply(brca.ciber.web[id.tumor,],2,function(x) sum(x>0)), MIXTURE=apply(GetMixture(brca.mix.n)[id.tumor,],2,function(x) sum(x>0)))/nrow(brca.ciber.web[id.tumor,]),2), file ="TablaOnlyTumor.txt")
#                               QUANTISEQ CIBERSORT MIXTURE
# B.cells.naive                    86.21     92.33   37.35
# B.cells.memory                    9.41     18.54    6.03
# Plasma.cells                     92.88     82.01   58.90
# T.cells.CD8                       8.86     72.15   32.88
# T.cells.CD4.naive                 0.18      1.83    0.27
# T.cells.CD4.memory.resting       68.04     95.98   37.81
# T.cells.CD4.memory.activated      1.28     19.00    3.56
# T.cells.follicular.helper        81.92     99.54   93.79
# T.cells.regulatory..Tregs.       82.92     88.49   48.68
# T.cells.gamma.delta               1.00      7.03    1.00
# NK.cells.resting                 34.61     78.26   12.97
# NK.cells.activated                0.37     41.10    3.01
# Monocytes                        53.15     53.97    8.40
# Macrophages.M0                   87.67     87.49   79.00
# Macrophages.M1                   94.25     92.42   77.63
# Macrophages.M2                   92.97     99.18   97.72
# Dendritic.cells.resting          78.36     27.58    8.22
# Dendritic.cells.activated        45.48     34.61   17.53
# Mast.cells.resting               89.86     96.89   49.95
# Mast.cells.activated              4.20      3.74    0.46
# Eosinophils                      46.76      5.39    0.00
# Neutrophils                      15.89     12.69    1.19

round(100*cbind(QUANTISEQ = apply(brca.quanti[id.tumor,-c(1,24)],2, function(x) sum(x>0,na.rm=T)),
                CIBERSORT=apply(brca.ciber.web[id.tumor,],2,function(x) sum(x>0)), 
                MIXTURE=apply(GetMixture(brca.mix.n)[id.tumor,],2,function(x) sum(x>0)))/nrow(brca.ciber.web[id.tumor,]),2)

M.mix <- GetMixture(brca.mix.n)[id.tumor,]
M.cib <- brca.ciber.web[id.tumor,]
M.q <- data.matrix(brca.quanti[id.tumor,-c(1,24)])
ct.names <- c("Bn","Bm","Pcs","CD8T","CD4T","CD4Tr","CD4Ta","Tfh","Tregs","gdt","Nkr","Nka","MM","M0","M1","M2","Dcr","Dca","Mcr","Mca","E","N")
df.bxp <- data.frame(betahat = c(as.numeric(M.mix), as.numeric(M.cib), as.numeric(M.q)),
                     CT = rep(ct.names, each = nrow(M.mix), times = 3),
                     Method = factor(rep(c("MIXTURE","CIBERSORT","QUANTISEQ"),each=prod(dim(M.mix))),
                                     levels = c("QUANTISEQ","CIBERSORT", "MIXTURE")))

LM22.bx <- ggplot(df.bxp, aes(x=CT, y=betahat, colour = Method)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs( x = " ",y = "Estimated Proportion")+ ggtitle("LM22") +facet_grid(~Method)                     


TIL10.M.mix <- GetMixture(TIL10.brca.mix.n)[id.tumor,]
TIL10.M.cib <- TIL10.brca.ciber.web[id.tumor,]
TIL10.M.q <- data.matrix(TIL10.brca.quanti[id.tumor,-c(1,12)])
colnames(TIL10.M.q)
TIL10.ct.names <- c("B","M1","M2","Mo","N","NK","CD4T","CD8T","Tregs","D")
TIL10.df.bxp <- data.frame(betahat = c(as.numeric(TIL10.M.mix), as.numeric(TIL10.M.cib), as.numeric(TIL10.M.q)),
                     CT = rep(TIL10.ct.names, each = nrow(TIL10.M.mix), times = 3),
                     Method = factor(rep(c("MIXTURE","CIBERSORT","QUANTISEQ"),each=prod(dim(TIL10.M.mix))),
                                     levels = c("QUANTISEQ","CIBERSORT", "MIXTURE")))

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
colnames(M.q)
ct.names <- c("Bn","Bm","Pcs","CD8T","CD4T","CD4Tr","CD4Ta","Tfh","Tregs","gdt","Nkr","Nka","MM","M0","M1","M2","Dcr","Dca","Mcr","Mca","E","N")

df.bxp <- data.frame(betahat = c(as.numeric(M.mix), as.numeric(M.cib), as.numeric(M.q)),
                     CT = rep(colnames(M.q), each = nrow(M.mix), times = 3),
                     Method = factor(rep(c("MIXTURE","CIBERSORT","QUANTISEQ"),each=prod(dim(M.mix))),
                                     levels = c("QUANTISEQ","CIBERSORT", "MIXTURE")))

LM22.bx <- ggplot(df.bxp, aes(x=CT, y=betahat, colour = Method)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+facet_grid(~Method)                     
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
CC.mat <- do.call(rbind, lapply(colnames(M.mix), function(x) {
  c(cor(M.mix[,x],TIL10.M.mix[,x], use = "pairwise.complete.obs"),
    cor(M.cib[,x],TIL10.M.cib[,x], use = "pairwise.complete.obs"),
    cor(M.q[,x],TIL10.M.q[,x], use = "pairwise.complete.obs"))
}))
rownames(CC.mat) <- colnames(M.mix)
colnames(CC.mat) <- c("MIXTURE","CIBERSORT","QUANTISEQ")
summary(CC.mat)
#           MIXTURE          CIBERSORT          QUANTISEQ          
# Min.   :-0.0002549   Min.   :0.0295   Min.   :-0.19288  
# 1st Qu.: 0.1130936   1st Qu.:0.1311   1st Qu.:-0.04651  
# Median : 0.2513331   Median :0.2231   Median :-0.01046  
# Mean   : 0.2608010   Mean   :0.2442   Mean   : 0.04274  
# 3rd Qu.: 0.3917043   3rd Qu.:0.3312   3rd Qu.: 0.17152  
# Max.   : 0.5135764   Max.   :0.5425   Max.   : 0.35687
pdf("/home/elmer/Dropbox/IDEAS/cibersort/GenomeR/Figuras/FigurasFinales/CorrMat_BRCA.pdf", paper = "a4", width = 5, height = 5)
Heatmap(CC.mat, cluster_rows = FALSE, cluster_columns = FALSE, name = "Correlation")
dev.off()
png("/home/elmer/Dropbox/IDEAS/cibersort/GenomeR/Figuras/FigurasFinales/CorrMat_BRCA.png")#, paper = "a4", width = 5, height = 5)
Heatmap(CC.mat, cluster_rows = FALSE, cluster_columns = FALSE, name = "Correlation")
dev.off()

round(CC.mat,3)
#         MIXTURE CIBERSORT QUANTISEQ
# B       0.482     0.390     0.357
# CD8T    0.514     0.543     0.223
# CD4T    0.219     0.178    -0.031
# Tregs   0.284     0.216    -0.052
# NK      0.072     0.115    -0.193
# Mo      0.091     0.308    -0.028
# M1      0.179     0.231     0.199
# M2      0.367     0.339    -0.145
# D       0.000     0.029     0.007
# N       0.400     0.094     0.091

wilcox.test(CC.mat[,"MIXTURE"],CC.mat[,"QUANTISEQ"], paired = T, alternative = "greater")
wilcox.test(CC.mat[,"CIBERSORT"],CC.mat[,"QUANTISEQ"], paired = T, alternative = "greater")
wilcox.test(CC.mat[,"MIXTURE"], CC.mat[,"CIBERSORT"], paired = T,  alternative = "greater")

wilcox.test(CC.mat[,"MIXTURE"],CC.mat[,"CIBERSORT"], paired = T, alternative = "g")


wilcox.test(CC.mat[c(1,2,3,4,8),"MIXTURE"],CC.mat[c(1,2,3,4,8),"CIBERSORT"], paired = T, alternative = "g")
# Wilcoxon signed rank test
# 
# data:  CC.mat[c(1, 2, 3, 4, 8), "MIXTURE"] and CC.mat[c(1, 2, 3, 4, 8), "CIBERSORT"]
# V = 13, p-value = 0.09375
# alternative hypothesis: true location shift is greater than 0
##Mormal samples
id.mormal <- which(brca$targets$sample =="11")


round(100*cbind(QUANTISEQ = apply(brca.quanti[id.mormal,-c(1,24)],2, function(x) sum(x>0,na.rm=T)),
  CIBERSORT=apply(brca.ciber.web[id.mormal,],2,function(x) sum(x>0)), MIXTURE=apply(GetMixture(brca.mix.n)[id.mormal,],2,function(x) sum(x>0)))/nrow(brca.ciber.web[id.mormal,]),2)
#                               CIBERSORT MIXTURE
# B cells naive                    87.61   11.50
# B cells memory                    4.42    1.77
# Plasma cells                     93.81   82.30
# T cells CD8                      71.68   22.12
# T cells CD4 naive                 0.88    0.00
# T cells CD4 memory resting      100.00   85.84
# T cells CD4 memory activated      1.77    0.00
# T cells follicular helper        94.69   63.72
# T cells regulatory  Tregs        38.05    0.88
# T cells gamma delta               2.65    0.00
# NK cells resting                 69.03   25.66
# NK cells activated               68.14   10.62
# Monocytes                        92.92   30.09
# Macrophages M0                   44.25   24.78
# Macrophages M1                   80.53   48.67
# Macrophages M2                  100.00   98.23
# Dendritic cells resting          28.32    1.77
# Dendritic cells activated        50.44   21.24
# Mast cells resting               94.69   57.52
# Mast cells activated              7.96    4.42
# Eosinophils                      18.58    0.00
# Neutrophils                       7.08    0.88


## ALl methods
CT.quanti <- data.matrix(brca.quanti[,-c(1,24)])
cellt <- "T cells follicular helper"
cellt <- "Macrophages M1"
acronim <- "M2"
clasif <- as.factor(brca$targets$pam50.pbcmc)
colnames(CT.quanti) <- colnames(GetMixture(brca.mix.n,"prop"))
cutoff <- rep(NA,length(colnames(GetMixture(brca.mix.n,"prop"))))
pval <- cutoff

colnames(CT.quanti) <- colnames(GetMixture(brca.mix.n,"prop"))
colnames(brca.ciber.web) <- colnames(GetMixture(brca.mix.n,"prop"))

#with LM22
Method.test <- "MIXTURE"
All <- lapply(colnames(CT.quanti), function(ct, sbt) {
  print(ct)
   
  df.m2 <- data.frame(Clasif = rep(clasif[id.tumor], 3),
                      CT = c(GetMixture(brca.mix.n,"prop")[id.tumor,ct],
                             brca.ciber.web[id.tumor,ct],
                             CT.quanti[id.tumor,ct]),
                      Method = factor(rep(c("MIXTURE", "CIBERSORT","QUANTISEQ"), each = length(id.tumor)),
                                      levels = c("QUANTISEQ","CIBERSORT","MIXTURE")),
                      Time = rep(brca$targets$survival.days[id.tumor] /356, 3),
                      Event = rep(brca$targets$vital.status[id.tumor], 3),
                      PBCMC = rep(brca$targets$pam50.pbcmc[id.tumor],3),
                      PAM50 = rep(brca$targets$pam50[id.tumor],3),
                      ER = rep(brca$targets$er[id.tumor],3) )
  if(is.null(sbt) ) {
    dm <- subset(df.m2, Method == Method.test )
  }else{
    dm <- subset(df.m2, Method == Method.test & PAM50 == sbt)
  }
  
  return(FindSurvCutoff(dm, 20, plot=F, title =colnames(GetMixture(brca.mix.n,"prop"))[ct]))
}, sbt = NULL)
names(All) <- colnames(CT.quanti)

sbts.PAM50 <- lapply(unique(brca$targets$pam50),function(sbt) {
  retL <- lapply(colnames(CT.quanti), function(ct, sbt) {
  # print(ct)
  df.m2 <- data.frame(Clasif = rep(clasif[id.tumor], 3),
                      CT = c(GetMixture(brca.mix.n,"prop")[id.tumor,ct],
                             brca.ciber.web[id.tumor,ct],
                             CT.quanti[id.tumor,ct]),
                      Method = factor(rep(c("MIXTURE", "CIBERSORT","QUANTISEQ"), each = length(id.tumor)),
                                      levels = c("QUANTISEQ","CIBERSORT","MIXTURE")),
                      Time = rep(brca$targets$survival.days[id.tumor] /356, 3),
                      Event = rep(brca$targets$vital.status[id.tumor], 3),
                      PBCMC = rep(brca$targets$pam50.pbcmc[id.tumor],3),
                      PAM50 = rep(brca$targets$pam50[id.tumor],3),
                      ER = rep(brca$targets$er[id.tumor],3) )
  if(is.null(sbt) ) {
    dm <- subset(df.m2, Method == Method.test )
  }else{
    print(sbt)
    dm <- subset(df.m2, Method == Method.test & PAM50 == sbt)
  }
  return(FindSurvCutoff(dm, 20, plot=F, title =colnames(GetMixture(brca.mix.n,"prop"))[ct]))
}, sbt = sbt)
 names(retL) <- colnames(CT.quanti)
 return(retL)
})
names(sbts.PAM50) <- unique(brca$targets$pam50)

sbts.PBCMC <- lapply(unique(brca$targets$pam50.pbcmc),function(sbt) {
  retL <- lapply(colnames(CT.quanti), function(ct, sbt) {
    # print(ct)
    df.m2 <- data.frame(Clasif = rep(clasif[id.tumor], 3),
                        CT = c(GetMixture(brca.mix.n,"prop")[id.tumor,ct],
                               brca.ciber.web[id.tumor,ct],
                               CT.quanti[id.tumor,ct]),
                        Method = factor(rep(c("MIXTURE", "CIBERSORT","QUANTISEQ"), each = length(id.tumor)),
                                        levels = c("QUANTISEQ","CIBERSORT","MIXTURE")),
                        Time = rep(brca$targets$survival.days[id.tumor] /356, 3),
                        Event = rep(brca$targets$vital.status[id.tumor], 3),
                        PBCMC = rep(brca$targets$pam50.pbcmc[id.tumor],3),
                        PAM50 = rep(brca$targets$pam50[id.tumor],3),
                        ER = rep(brca$targets$er[id.tumor],3) )
    if(is.null(sbt) ) {
      dm <- subset(df.m2, Method == Method.test )
    }else{
      dm <- subset(df.m2, Method == Method.test & PBCMC == sbt)
    }
    return(FindSurvCutoff(dm, 20, plot=F, title =colnames(GetMixture(brca.mix.n,"prop"))[ct]))
  }, sbt = sbt)
  names(retL) <- colnames(CT.quanti)
  return(retL)
})
names(sbts.PBCMC) <- unique(brca$targets$pam50.pbcmc)

sbts.ER <- lapply(unique(brca$targets$er)[1:2],function(sbt) {
  retL <- lapply(colnames(CT.quanti), function(ct, sbt) {
    # print(ct)
    df.m2 <- data.frame(Clasif = rep(clasif[id.tumor], 3),
                        CT = c(GetMixture(brca.mix.n,"prop")[id.tumor,ct],
                               brca.ciber.web[id.tumor,ct],
                               CT.quanti[id.tumor,ct]),
                        Method = factor(rep(c("MIXTURE", "CIBERSORT","QUANTISEQ"), each = length(id.tumor)),
                                        levels = c("QUANTISEQ","CIBERSORT","MIXTURE")),
                        Time = rep(brca$targets$survival.days[id.tumor] /356, 3),
                        Event = rep(brca$targets$vital.status[id.tumor], 3),
                        PBCMC = rep(brca$targets$pam50.pbcmc[id.tumor],3),
                        PAM50 = rep(brca$targets$pam50[id.tumor],3),
                        ER = rep(brca$targets$er[id.tumor],3) )
    if(is.null(sbt) ) {
      dm <- subset(df.m2, Method == Method.test )
    }else{
      dm <- subset(df.m2, Method == Method.test & ER == sbt)
    }
    return(FindSurvCutoff(dm, 20, plot=F, title =colnames(GetMixture(brca.mix.n,"prop"))[ct]))
  }, sbt = sbt)
  names(retL) <- colnames(CT.quanti)
  return(retL)
})
names(sbts.ER) <- unique(brca$targets$er)[1:2]

# return(list(gplot=gs, cutoff = ic, pval = pval, ssum <- surv_summary(fit.m)))
pval <- unlist(lapply(All, function(x) x$pval))
All.ct.surv <- rep(NA, 22)
All.ct.surv[pval < 0.05] <- unlist(lapply(names(All)[pval < 0.05], function(x) if(diff(attr(All[[x]]$surv.summary,"table")$median) > 0){
  "L"
}else "H"))

PAM50.ct.surv <- do.call(cbind, lapply(sbts.PAM50, function(s){
  vv <- rep(NA,22)
  names(vv) <- colnames(CT.quanti)
  val <- unlist(lapply(s, function(x) x$pval))
  vv[names(vv) %in% names(val)] <- val
  vv.sum <- rep(NA,22)
  vv.sum[which(vv < 0.05)] <- unlist(lapply(names(s)[which(vv < 0.05)], function(x) if(diff(attr(s[[x]]$surv.summary,"table")[,"*rmean"]) > 0){
    "L"
  }else "H"))
  return(vv.sum)
}
  ))

PBCMC.ct.surv <- do.call(cbind, lapply(sbts.PBCMC, function(s){
  
  vv <- rep(NA,22)
  names(vv) <- colnames(CT.quanti)
  val <- unlist(lapply(s, function(x) x$pval))
  
  vv[names(vv) %in% names(val)] <- val
  vv.sum <- rep(NA,22)
  vv.sum[which(vv < 0.05)] <- unlist(lapply(names(s)[which(vv < 0.05)], function(x) if(diff(attr(s[[x]]$surv.summary,"table")[,"*rmean"]) > 0){
    "L"
  }else "H"))
  return(vv.sum)
}
))
ER.ct.surv <- do.call(cbind, lapply(sbts.ER, function(s){
  
  vv <- rep(NA,22)
  names(vv) <- colnames(CT.quanti)
  val <- unlist(lapply(s, function(x) x$pval))
  vv[names(vv) %in% names(val)] <- val
  vv.sum <- rep(NA,22)
  vv.sum[which(vv < 0.05)] <- unlist(lapply(names(s)[which(vv < 0.05)], function(x) if(diff(attr(s[[x]]$surv.summary,"table")[,"*rmean"]) > 0){
    "L"
  }else "H"))
  
  return(vv.sum)
}
))
CT.surv <- cbind(All.ct.surv, PAM50.ct.surv,ER.ct.surv)
CT.surv <- cbind(All.ct.surv, PBCMC.ct.surv,ER.ct.surv)
CT.surv <- CT.surv[,-5]##Normal left out
CT.surv <- CT.surv[,-4]##Normal left out
colnames(CT.surv) <- c("All", "Not Assigned","Her2","Luminal B","Luminal A","Basal","ER+","ER-")
colnames(CT.surv) <- c("All","Luminal B","Her2","Luminal A","Basal","ER+","ER-")

rownames(CT.surv) <- colnames(CT.quanti)
View(CT.surv)
pdf(paste("/home/elmer/Dropbox/IDEAS/cibersort/GenomeR/Figuras/Survival_BRCA",Method.test,".pdf",sep=""))
Heatmap(CT.surv, cluster_rows = FALSE, cluster_columns = FALSE, col = c("L"="green","H"="red"), 
        na_col ="white", rect_gp = gpar(col = "black", lwd = 1),column_names_side = "top",name ="Prognosis Association")
dev.off()


table(df.m2$PAM50)

# dm <- subset(df.m2, Method == "CIBERSORT")

TIL10.quanti <- TIL10.brca.quanti2[,-c(1,12)]


colnames(TIL10.brca.ciber.web) <- colnames(TIL10.quanti)
TIL10.cellt <- "Macrophages.M2"
TIL10.df.m2 <- data.frame(Clasif = rep(clasif[id.tumor], 3),
                    CT = c(GetMixture(TIL10.brca.mix.n,"prop")[id.tumor,TIL10.cellt],
                           TIL10.brca.ciber.web[id.tumor,TIL10.cellt],
                           TIL10.quanti[id.tumor,TIL10.cellt]),
                    Method = factor(rep(c("MIXTURE", "CIBERSORT","QUANTISEQ"), each = length(id.tumor)),
                                    levels = c("QUANTISEQ","CIBERSORT","MIXTURE")),
                    Time = rep(brca$targets$survival.days[id.tumor] /356, 3),
                    Event = rep(brca$targets$vital.status[id.tumor], 3),
                    PBCMC = rep(brca$targets$pam50.pbcmc[id.tumor],3),
                    PAM50 = rep(brca$targets$pam50[id.tumor],3),
                    ER = rep(brca$targets$er[id.tumor],3)
)


for(s in c("Basal","Her2","LumB","LumA")){
  lapply( c("CIBERSORT","MIXTURE"), function(i){
    df.test.mix <- subset(df.m2, Method == i & PBCMC == s   )
    df.test.mix$Q <- quantilize(df.test.mix$CT, c(0.33,0.66))
    fit.m <- survfit(Surv(Time, Event) ~ Q, data = df.test.mix)
    gs <- ggsurvplot(fit.m, df.test.mix, pval = TRUE) + ggtitle(paste(i,s,cellt))
    print(table(df.test.mix$Q, df.test.mix$PBCMC))
    #return(gs)
    print(gs)
  })  
}


gs.tfh <- lapply( c("QUANTISEQ","CIBERSORT","MIXTURE"), function(i){
  df.test.mix <- subset(df.m2, Method == i )
  df.test.mix$Q <- quantilize(df.test.mix$CT, c(0.33,0.66))
  fit.m <- survfit(Surv(Time, Event) ~ Q, data = df.test.mix)
  gs <- ggsurvplot(fit.m, df.test.mix, pval = TRUE, 
                   legend.labs = paste(acronim,c("Low", "Midle", "High")), xlab = "Time (months)") + ggtitle(paste(i,"(",cellt,")"))
  print(table(df.test.mix$Q, df.test.mix$PBCMC))
  return(gs)
  #print(gs)
})



gs.M2 <- lapply( c("QUANTISEQ","CIBERSORT","MIXTURE"), function(i){
  df.test.mix <- subset(df.m2, Method == i )
   df.test.mix$Q <- quantilize(df.test.mix$CT, c(0.33,0.66))
   df.test.mix$Q[df.test.mix$Q == "Q2"] <- "Q3"
   df.test.mix$Q <- factor(as.character(df.test.mix$Q))
  fit.m <- survfit(Surv(Time, Event) ~ Q, data = df.test.mix)
  gs <- ggsurvplot(fit.m, df.test.mix, pval = TRUE)#, 
                   # legend.labs = paste(acronim,c("Low", "Midle", "High")), xlab = "Time (months)") + ggtitle(paste(i,"(",cellt,")"))
  print(table(df.test.mix$Q, df.test.mix$PBCMC))
  return(gs)
  #print(gs)
})

TIL10.gs.M2 <- lapply( c("QUANTISEQ","CIBERSORT","MIXTURE"), function(i){
  df.test.mix <- subset(TIL10.df.m2, Method == i )
  df.test.mix$Q <- quantilize(df.test.mix$CT, c(0.33,0.66))
  df.test.mix$Q <- quantilize(df.test.mix$CT, c(0.6))
  df.test.mix$Q[df.test.mix$Q == "Q2"] <- "Q3"
  df.test.mix$Q <- factor(as.character(df.test.mix$Q))
  fit.m <- survfit(Surv(Time, Event) ~ Q, data = df.test.mix)
  gs <- ggsurvplot(fit.m, df.test.mix, pval = TRUE)#, 
                   # legend.labs = paste(acronim,c( "Low","Midle", "High")), xlab = "Time (months)") + 
    ggtitle(paste(i,"(",TIL10.cellt,")"))
  print(table(df.test.mix$Q, df.test.mix$PBCMC))
  return(gs)
  #print(gs)
})
# setEPS()
# postscript(file = "BRCA_Full_M1_Surv.eps")
 arrange_ggsurvplots(gs.tfh, print = TRUE,
                     ncol = 3, nrow = 1, risk.table.height = 0.4)
 
# dev.off()
 arrange_ggsurvplots(gs.M2, print = TRUE,
                     ncol = 3, nrow = 1, risk.table.height = 0.4)
 arrange_ggsurvplots(TIL10.gs.M2, print = TRUE,
                     ncol = 3, nrow = 1, risk.table.height = 0.4)
 gs.aux <- gs.tfh
 
length(gs.aux)
gs.aux[[4]] <- TIL10.gs.M2[[1]]
gs.aux[[5]] <- TIL10.gs.M2[[2]]
gs.aux[[6]] <- TIL10.gs.M2[[3]] 
class(gs.aux)                 
 length(gs.aux)                 
 arrange_ggsurvplots(gs.aux, print = TRUE,
                     ncol = 3, nrow =2, risk.table.height = 0.4)
 
 # 
# 
# pdf(file = paste("Survival_BRCA_",acronim,cellt,".pdf",sep=""), width = 14, paper = "a4r")
# arrange_ggsurvplots(gs, print = TRUE,
#                     ncol = 2, nrow = 1, risk.table.height = 0.4)
# dev.off()
# 
# jpeg(file = paste("Survival_BRCA_",acronim,cellt,".jpg",sep=""), width = 600)
# arrange_ggsurvplots(gs, print = TRUE,
#                     ncol = 2, nrow = 1, risk.table.height = 0.4)
# dev.off()
# 

CT.mix.presence <- GetCellTypes(brca.mix.n)


for(i in c("CIBERSORT","MIXTURE")){
  df.test.mix <- subset(df.m2, Method == i  & ER == "negative" )
  df.test.mix$qe <- quantilize(df.test.mix$CT, c(0.33,0.66))
  fit.m <- survfit(Surv(Time, Event) ~ qe, data = df.test.mix)
  gs <- ggsurvplot(fit.m, df.test.mix, pval = TRUE) + ggtitle(paste(i, "ER-", cellt))
  print(gs)
}

lapply( c("CIBERSORT","MIXTURE"), function(i){
  df.test.mix <- subset(df.m2, Method == i  & ER == "positive" )
  df.test.mix$qe <- quantilize(df.test.mix$CT,c(0.33,0.66))
  fit.m <- survfit(Surv(Time, Event) ~ qe, data = df.test.mix)
  gs <- ggsurvplot(fit.m, df.test.mix, pval = TRUE) + ggtitle(paste(i, "ER+",cellt))
  print(gs)
}
)


table(df.test.mix$qe)
table(df.m2$PBCMC)
quantile(df.test.mix$CT, c(0.25,0.5,0.75))

df.test.cib <- subset(df.m2, Method %in% c("CIBERSORT") )
df.test.cib$qe <- quartilize(df.test.cib$CT)

# df.test.cib$qe[df.test.cib$CT > quantile(df.test.cib$CT, 0.33)] <- "M-CIB"
# df.test.cib$qe[df.test.cib$CT > quantile(df.test.cib$CT, 0.66)] <- "H-CIB"
fit <- survfit(Surv(Time, Event) ~ qe, data = df.test.cib)
ggsurvplot(fit, df.test.cib, pval = TRUE)+ ggtitle("CIBERSORT")

for(i in levels(df.m2$Method )){
  df.test.mix <- subset(df.m2, Method == i   )
  df.test.mix$qe <- quartilize(df.test.mix$CT)
  # df.test.mix$qe <- "L-MIX"
  # df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.33)] <- "M-MIX"
  # df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.66)] <- "H-MIX"
  fit.m <- survfit(Surv(Time, Event) ~ qe, data = df.test.mix)
  gs <- ggsurvplot(fit.m, df.test.mix, pval = TRUE) + ggtitle(i)
  print(gs)
}


for(i in c("CIBERSORT","MIXTURE")){
  df.test.mix <- subset(df.m2, Method == i   )
  df.test.mix$qe <- quartilize(df.test.mix$CT)
  # df.test.mix$qe <- "L-MIX"
  # df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.33)] <- "M-MIX"
  # df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.66)] <- "H-MIX"
  fit.m <- survfit(Surv(Time, Event) ~ qe, data = df.test.mix)
  gs <- ggsurvplot(fit.m, df.test.mix, pval = TRUE) + ggtitle(i)
  print(gs)
}

df.test.mix <- subset(df.m2, Method %in% c("MIXTURE")   )
df.test.mix$qe <- quartilize(df.test.mix$CT)
# df.test.mix$qe <- "L-MIX"
# df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.33)] <- "M-MIX"
# df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.66)] <- "H-MIX"
fit.m <- survfit(Surv(Time, Event) ~ qe, data = df.test.mix)
ggsurvplot(fit.m, df.test.mix, pval = TRUE) + ggtitle("MIXTURE")


df.test.cib <- subset(df.m2, Method %in% c("CIBERSORT") & PBCMC == "Basal" )
df.test.cib$qe <- "L-CIB"
df.test.cib$qe[df.test.cib$CT > quantile(df.test.cib$CT, 0.33)] <- "M-CIB"
df.test.cib$qe[df.test.cib$CT > quantile(df.test.cib$CT, 0.66)] <- "H-CIB"
fit <- survfit(Surv(Time, Event) ~ qe, data = df.test.cib)
ggsurvplot(fit, df.test.cib, pval = TRUE)

df.test.mix <- subset(df.m2, Method %in% c("MIXTURE")   & PBCMC == "Basal")
df.test.mix$qe <- "L-MIX"
df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.33)] <- "M-MIX"
df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.66)] <- "H-MIX"
fit <- survfit(Surv(Time, Event) ~ qe, data = df.test.mix)
ggsurvplot(fit, df.test.mix, pval = TRUE)



IType <- "Her2"
df.test.cib <- subset(df.m2, Method %in% c("CIBERSORT") & PBCMC == IType )
df.test.cib$qe <- "L-CIB"
df.test.cib$qe[df.test.cib$CT > quantile(df.test.cib$CT, 0.33)] <- "M-CIB"
df.test.cib$qe[df.test.cib$CT > quantile(df.test.cib$CT, 0.66)] <- "H-CIB"
fit <- survfit(Surv(Time, Event) ~ qe, data = df.test.cib)
ggsurvplot(fit, df.test.cib, pval = TRUE)

df.test.mix <- subset(df.m2, Method %in% c("MIXTURE")   & PBCMC == IType)
df.test.mix$qe <- "L-MIX"
df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.33)] <- "M-MIX"
df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.66)] <- "H-MIX"
fit <- survfit(Surv(Time, Event) ~ qe, data = df.test.mix)
ggsurvplot(fit, df.test.mix, pval = TRUE)


IType <- "LumB"
df.test.cib <- subset(df.m2, Method %in% c("CIBERSORT") & PBCMC == IType )
df.test.cib$qe <- "L-CIB"
df.test.cib$qe[df.test.cib$CT > quantile(df.test.cib$CT, 0.33)] <- "M-CIB"
df.test.cib$qe[df.test.cib$CT > quantile(df.test.cib$CT, 0.66)] <- "H-CIB"
fit <- survfit(Surv(Time, Event) ~ qe, data = df.test.cib)
ggsurvplot(fit, df.test.cib, pval = TRUE)

df.test.mix <- subset(df.m2, Method %in% c("MIXTURE")   & PBCMC == IType)
df.test.mix$qe <- "L-MIX"
df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.33)] <- "M-MIX"
df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.66)] <- "H-MIX"
fit <- survfit(Surv(Time, Event) ~ qe, data = df.test.mix)
ggsurvplot(fit, df.test.mix, pval = TRUE)

IType <- "LumA"
df.test.cib <- subset(df.m2, Method %in% c("CIBERSORT") & PBCMC == IType )
df.test.cib$qe <- "L-CIB"
df.test.cib$qe[df.test.cib$CT > quantile(df.test.cib$CT, 0.33)] <- "M-CIB"
df.test.cib$qe[df.test.cib$CT > quantile(df.test.cib$CT, 0.66)] <- "H-CIB"
fit <- survfit(Surv(Time, Event) ~ qe, data = df.test.cib)
ggsurvplot(fit, df.test.cib, pval = TRUE)

df.test.mix <- subset(df.m2, Method %in% c("MIXTURE")   & PBCMC == IType)
df.test.mix$qe <- "L-MIX"
df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.33)] <- "M-MIX"
df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.66)] <- "H-MIX"
fit <- survfit(Surv(Time, Event) ~ qe, data = df.test.mix)
ggsurvplot(fit, df.test.mix, pval = TRUE)



for(i in unique(brca$targets$pam50.pbcmc)){
  df.test.mix <- subset(df.m2, Method %in% c("MIXTURE")   & PBCMC == i)
  df.test.mix$qe <- "L-MIX"
  df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.33)] <- NA
  df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.66)] <- "H-MIX"
  fit <- survfit(Surv(Time, Event) ~ qe, data = df.test.mix)
  print(ggsurvplot(fit, df.test.mix, pval = TRUE))
  
}


for(i in levels(df.m2$Method)){
  df.test.mix <- subset(df.m2, Method == i)
  df.test.mix$qe <- "L"
  df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.33)] <- "M"
  df.test.mix$qe[df.test.mix$CT > quantile(df.test.mix$CT, 0.66)] <- "H"
  fit <- survfit(Surv(Time, Event) ~ qe, data = df.test.mix)
  print(ggsurvplot(fit, df.test.mix, pval = TRUE))
  
}

df.brca.cm <- rbind(df.test.cib, df.test.mix)


fit <- survfit(Surv(Time, Event) ~ qe, data = df.brca.cm)
ggsurvplot(fit, df.brca.cm, pval = TRUE)


fit <- survfit(Surv(Time, Event) ~ qe, data = subset(df.brca.cm, qe %in% c("L-CIB","L-MIX")))
ggsurvplot(fit, subset(df.brca.cm, qe %in% c("L-CIB","L-MIX")), pval = TRUE)



ggplot(na.omit(df.m2)) +
  geom_boxplot(aes(x = Clasif, y = CT, fill= Method))  + xlab("(survival")




wilcox.test(CT~Clasif, subset(df.m2, Method == "CIBERSORT"))
# Wilcoxon rank sum test with continuity correction
# 
# data:  CT by Clasif
# W = 8116, p-value = 0.0245
# alternative hypothesis: true location shift is not equal to 0  
wilcox.test(CT~Clasif, subset(df.m2, Method == "MIXTURE"))

# Wilcoxon rank sum test with continuity correction
# 
# data:  CT by Clasif
# W = 6984, p-value = 0.0001235
# alternative hypothesis: true location shift is not equal to 0
wilcox.test(CT~Clasif, subset(df.m2, Method == "ABBAS"))
# Wilcoxon rank sum test with continuity correction
# 
# data:  CT by Clasif
# W = 9144, p-value = 0.412
# alternative hypothesis: true location shift is not equal to 0




# saveRDS(list(brca.cib, brca.mix2,brca.abbas), "MIXTURE_ABBAS_Results.rds")
