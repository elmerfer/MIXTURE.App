#We test here the TCGA with MIXER and CIBERSORT.
#

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

source('/Utils/MIXTURE.DEBUG_V0.1.R')
##change the directory to your own directory!!!
#the BRCA RNAseq data can be downloaded from https://www.dropbox.com/s/zki1gkx5mq1quah/BRCA_rna.rds?dl=0
# brca <- readRDS("/home/elmer/Dropbox/Doctorandos/DarioRocha/BRCA/processed_data/BRCA_rna.rds")
TNBC <- apply(brca$targets[,c("er","pgr","her2")],1, FUN = function(x) all(x == "negative"))
brca$targets$TNBC <- TNBC ##defining Triple Negative BRCA

quantilize <- function(var, qs = c(0.25,0.5,0.75)){
  varq <- quantile(var,qs)
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
 dge <- calcNormFactors(dge)
 brca.norm <- brca
 brca.norm$E <- cpm(dge$counts)
# 
 M.brca.n <- brca.norm$E
 sum(rownames(M.brca.n) %in% rownames(LM22))
rownames(M.brca.n)[1:3]
#df.brca.to.save <- data.frame(genes = rownames(M.brca.n)[rownames(M.brca.n) %in% rownames(LM22)], M.brca.n[rownames(M.brca.n) %in% rownames(LM22),])
#colnames(df.brca.to.save)[1] <- "Gene symbol"
# write.table(df.brca.to.save, file = "/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/BRCA.TCGA.txt", quote = F, row.names = F, sep= "\t")
#write.xlsx(df.brca.to.save, file = "/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/BRCA.TCGA.xlsx")

ncores2use <- 10L
# 
 # brca.cib.n <- MIXTURE(expressionMatrix = M.brca.n, signatureMatrix =  LM22, functionMixture =  cibersort, useCores = ncores2use, verbose  =  TRUE,
 #                       iter = 1000, nullDist = "PopulationBased")
# 
 brca.mix.n <- MIXTURE(expressionMatrix = M.brca.n, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = ncores2use, verbose  =  TRUE,
                       iter = 1000, nullDist = "PopulationBased")
# 
# 
 brca.abbas.n <- MIXTURE(expressionMatrix = M.brca.n, signatureMatrix =  LM22, functionMixture = ls.rfe.abbas, useCores = ncores2use, verbose  =  TRUE,
                         iter = 1000, nullDist = "none")
 
 
 Y.dt <- t(cpm(dge$counts, log =TRUE)[rownames(M.brca.n) %in% rownames(LM22),])
 dim(Y.dt)
 R.dt <- t(log2(LM22)[rownames(LM22) %in% colnames(Y.dt),])
 dim(R.dt)
 brca.dt.n <- dtangle(Y =  Y.dt , reference = R.dt)
 head(brca.dt.n$estimates[,1:10])
# 
brca.ciber.web <- ReadCibersortWebResults("Data/BRCA_CIBERSORT.csv", type = "csv")

#ABIS
brca.abis <- MIXTURE(expressionMatrix = M.brca.n, signatureMatrix =  LM22, functionMixture =  rlm.abis, useCores = ncores2use, verbose  =  TRUE,
                     iter = 1000, nullDist = "none")


#no norm                      
df.brca <- data.frame( Nb = c(apply(GetCellTypes(brca.mix.n),1, sum),
                              apply(GetCellTypes(brca.abbas.n),1, sum),
                              apply(brca.dt.n$estimates>0,1,sum),
                              apply(brca.ciber.web>0,1,sum),
                              apply(GetCellTypes(brca.abis),1, sum)),
                       method = factor(rep(c("ABBAS","MIXTURE","DTANGLE","CIBERSORT", "ABIS"),each = ncol(M.brca.n)), 
                                       levels = c("ABBAS","ABIS", "DTANGLE","CIBERSORT", "MIXTURE")))

ddply(df.brca, .(method), summarise, min= min(Nb, na.rm=T), Q1 = quantile(Nb, 0.25),  Q2 = median(Nb, na.rm=T), Q3 = quantile(Nb, 0.75), max = max(Nb, na.rm=T))
#      method min Q1 Q2 Q3 max
# 1     ABBAS   2  5  7  8  12
# 2      ABIS   7 11 12 13  17
# 3   DTANGLE  22 22 22 22  22
# 4 CIBERSORT   8 11 12 13  17
# 5   MIXTURE   0  4  5  7  11

##Ploting the amount of estimated cell - types per subject for each method
ggplot(df.brca, aes(x=method, y=Nb, fill=method)) +
  geom_boxplot(position=position_dodge(1)) + ggtitle("BRCA")

##Mixture cell types
CT.brca.mix <- GetMixture(brca.mix.n,"proportion")
CT.brca.abbas <- GetMixture(brca.abbas.n,"proportion")
CT.brca.dt <- brca.dt.n$estimates
CT.brca.abis <- GetMixture(brca.abis,"proportion")


##Survival by PAM50 intrinsic sutype (defined by PBCMC - Fresno et al. https://doi.org/10.1093/bioinformatics/btw704)
clasif <- as.factor(brca$targets$pam50.pbcmc)
#Short vs Long survival on Macrophages M2

# cellt <- "Macrophages M2"
#cellt <- "Plasma cells"
# df.m2 <- data.frame(Clasif = rep(clasif, 5),
#                     CT = c(CT.brca.dt[,cellt],
#                            GetMixture(brca.mix.n,"prop")[,cellt],
#                            GetMixture(brca.abbas.n,"prop")[,cellt],
#                            brca.ciber.web[,cellt],
#                            GetMixture(brca.abis,"prop")[,cellt]),
#                  Method = factor(rep(c("DTANGLE","MIXTURE", "ABBAS", "CIBERSORT", "ABIS"), each = ncol(brca$E)),
#                                  levels = c("ABBAS","ABIS","DTANGLE","CIBERSORT","MIXTURE")),
#                  Time = rep(brca$targets$survival.days /356, 5),
#                  Event = rep(brca$targets$vital.status, 5),
#                  PBCMC = rep(brca$targets$pam50.pbcmc,5),
#                  PAM50 = rep(brca$targets$pam50,5),
#                  ER = rep(brca$targets$er,5)
#                             )
# colnames(GetMixture(brca.mix.n,"prop"))

## Only CIBERSORT and MIXTURE
cellt <- "T cells follicular helper"
cellt <- "Macrophages M0"
acronim <- "M0"
df.m2 <- data.frame(Clasif = rep(clasif, 2),
                    CT = c(GetMixture(brca.mix.n,"prop")[,cellt],
                           brca.ciber.web[,cellt]),
                    Method = factor(rep(c("MIXTURE", "CIBERSORT"), each = ncol(brca$E)),
                                    levels = c("CIBERSORT","MIXTURE")),
                    Time = rep(brca$targets$survival.days /356, 2),
                    Event = rep(brca$targets$vital.status, 2),
                    PBCMC = rep(brca$targets$pam50.pbcmc,2),
                    PAM50 = rep(brca$targets$pam50,2),
                    ER = rep(brca$targets$er,2)
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


gs <- lapply( c("CIBERSORT","MIXTURE"), function(i){
  df.test.mix <- subset(df.m2, Method == i )
  df.test.mix$Q <- quantilize(df.test.mix$CT, c(0.33,0.66))
  fit.m <- survfit(Surv(Time, Event) ~ Q, data = df.test.mix)
  gs <- ggsurvplot(fit.m, df.test.mix, pval = TRUE, 
                   legend.labs = paste(acronim,c("Low", "Midle", "High")), xlab = "Time (months)") + ggtitle(paste(i,"(",cellt,")"))
  print(table(df.test.mix$Q, df.test.mix$PBCMC))
  return(gs)
  #print(gs)
})
# setEPS()
# postscript(file = "BRCA_Full_M1_Surv.eps")
# arrange_ggsurvplots(gs, print = TRUE,
#                     ncol = 2, nrow = 1, risk.table.height = 0.4)
# dev.off()
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


round(100*cbind(CIBERSORT=apply(brca.ciber.web,2,function(x) sum(x>0)), MIXTURE=apply(GetMixture(brca.mix.n),2,function(x) sum(x>0)))/nrow(brca.ciber.web),2)
write.table(round(100*cbind(CIBERSORT=apply(brca.ciber.web,2,function(x) sum(x>0)), MIXTURE=apply(GetMixture(brca.mix.n),2,function(x) sum(x>0)))/nrow(brca.ciber.web),2), file ="Tabla.txt")
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
