##REAL DATA from Noewman
rm(list=ls())
library(dtangle)
library(immunedeconv)
source('Utils/MIXTURE.DEBUG_V0.1.R')
##since FL only contains B, CD8 and CD4 cells, we merge (sum) the proportions of such cell according to Newman et al.
cell.types.names11 <- c("B","B","PC","CD8","CD4","CD4","CD4","CD4","CD4","TGD","NK","NK","Mo","Ma","Ma","Ma","D","D","Mt","Mt","Eo","N")

#Representacion de los tiposcelulares de LM22 para el ejemplo de immunodeconv
# Type	No.ct	B cells naive	B cells memory	Plasma cells	T cells CD8	T cells CD4 naive	T cells CD4 memory resting	T cells CD4 memory activated	T cells follicular helper	T cells regulatory (Tregs)	T cells gamma delta	NK cells resting	NK cells activated	Monocytes	Macrophages M0	Macrophages M1	Macrophages M2	Dendritic cells resting	Dendritic cells activated	Mast cells resting	Mast cells activated	Eosinophils	Neutrophils
# LM22-original	22 = c("Bn","Bm","Pcs","CD8T","CD4T","CD4Tr","CD4Ta","Tfh","Tregs","gdt","Nkr","Nka","MM","M0","M1","M2","Dcr","Dca","Mcr","Mca","E","N")
# LM22-merge4	4	   = c("B", "B", "B",  "CD8T","CD4T","CD4T", "CD4T", "CD4T","CD4T","CD8T","R",  "R",  "R", "R", "R","R",  "R",  "R",  "R",  "R", "R","R")
# then ctClassType = c("B","B","B","CD8T","CD4T","CD4T","CD4T","CD4T","CD4T","CD8T","R","R","R","R","R","R","R","R","R","R","R","R")
# TIL10 example	   = c("B", "B", "B",  "CD8T","CD4T","CD4T", "CD4T", "CD4T","CD4T","CD8T","R",  "R",  "R", "R", "R","R",  "D",  "D",  "R",  "R", "R","R")
#this represent a vector of 4 innmuno cell types group, as in https://www.nature.com/articles/s41587-019-0114-2
# LM22-original	22 = c("Bn","Bm","Pcs","CD8T","CD4T","CD4Tr","CD4Ta","Tfh","Tregs","gdt","Nkr","Nka","MM","M0","M1","M2","Dcr","Dca","Mcr","Mca","E","N")
LM22.til10.example <- c("R", "R", "R",  "CD8T","CD4T","CD4T", "CD4T", "CD4T","CD4T","CD8T","NK","NK","R","R","R","R", "D",  "D",  "R",  "R", "R","R")
##Follicular Data set only have fl.composite <- c("CD4T", "CD8T","B","R","R","R","R","R","R","R","R","R")
# LM22-original	22 = c("Bn","Bm","Pcs","CD8T","CD4T","CD4Tr","CD4Ta","Tfh","Tregs","gdt","Nkr","Nka","MM","M0","M1","M2","Dcr","Dca","Mcr","Mca","E","N")
LM22.to.FLtypes <-    c("B cells", "B cells", "B cells",  "CD8T","CD4T","CD4T", "CD4T", "CD4T","CD4T","CD8T","R","R","R","R","R","R", "R",  "R",  "R",  "R", "R","R") 
LM22.to.PBMCtypes <-  c("B cells naive", "B cells memory", "R",  "T cells CD8","T cells CD4 naive","R", "T cells CD4 memory activated", "R",
                        "R","T cells gamma delta","R","NK cells activated","Monocytes","R","R","R", "R",  "R",  "R",  "R", "R","R") 
  # [1] "B cells naive"                "B cells memory"               "T cells CD8"                  "T cells CD4 naive"            "T cells CD4 memory resting"  
# [6] "T cells CD4 memory activated" "T cells gamma delta"          "NK cells activated"           "Monocytes" 
#"B" "CD8" "CD4" "NK"  "Mo"  "M0"  "M1"  "M2"  "D"   "M"   "E"   "N" c
cbind(colnames(LM22), LM22.to.FLtypes)
#Representacion de los tipos celulares para TIL10
# colnames(TIL10)
# [1] "B.cells"         "Macrophages.M1"  "Macrophages.M2"  "Monocytes"       "Neutrophils"     "NK.cells"        "T.cells.CD4"    
# [8] "T.cells.CD8"     "Tregs"           "Dendritic.cells"
TIL10.ct <- c("B","M1","M2","Mo","N","NK","CD4T","CD8T","Tregs","D")
TIL10.ct.2 <- c("R","R","R","R","R","NK","CD4T","CD8T","CD4T","D")

CompositeMixture <- function(mixt, composite){
  colnames(mixt) <- composite
  
  ##This function merges and summarize the cell types proportion from the 22 cell-types to the 11 ones defined by Newman et al.  
  # dat <- data.matrix(t(aggregate(data.frame(t(mixt)), by= list(composite), FUN = sum)[,-1]))
  ret<- do.call(cbind, lapply(unique(composite), function(x) rowSums(mixt[,composite==x,drop=FALSE],na.rm = T)))
  colnames(ret) <- unique(composite)
  return(ret)
  
}



##load Data from Newman (downloaded from DTANGLE source code)
load("Data/newman_pbmc.rda")
load("Data/newman_fl.rda")
load("Data/LM22.RData")
##prepare data for CIBERSORT and MIXTURE
FL <- t(2^newman_fl$data$log)[,1:14]
PBMC <- t(2^newman_pbmc$data$log)[,1:20]

##write data for CIBERSORT web site
# write.table(FL, file="/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/NewmanFL.txt", quote = FALSE, row.names = TRUE, sep="\t")
# write.table(PBMC, file="/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/NewmanPBCMC.txt", quote = FALSE, row.names = TRUE, sep="\t")

rn <- rownames(PBMC)

# for( i in 1:20){
#   df.out <- data.frame("Gene symbol" = rn)
#   df.out$pp <- PBMC[,i]
#   colnames(df.out) <- c("Gene symbol", paste("PBMC_",i,sep=""))
#   write.table(df.out, file=paste("/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/",colnames(df.out)[2],".csv",sep=""), quote = FALSE, row.names = FALSE, sep="\t")
# }


##Read Results from CIBERSORT web site
FL.cib.site <- ReadCibersortWebResults(file = "Data/FL_CIBERSORT.Output_Job6.csv")
PBMC.cib.site <- ReadCibersortWebResults(file = "Data/PBMC_CIBERSORT.Output_Job7.csv")

##read single subject CIBERSORT web

# FL.cib.ss <- do.call(rbind, lapply(1:14, function(x) {
#   print(paste("FL",x))
#   cat("\n")
#   read.csv(paste("/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/CIBres/FL",x,".csv", sep=""),h=T)
# } ))
# colnames(FL.cib.ss) <- str_replace_all(colnames(FL.cib.ss),"\\.", " " )
# colnames(FL.cib.ss)[10] <- "T cells regulatory  (Tregs)" 
# rownames(FL.cib.ss) <- paste("FL",1:14,sep=" ")
# 
# PBMC.cib.ss <- do.call(rbind, lapply(1:20, function(x) {
#   print(paste("PBMC",x))
#   read.csv(paste("/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/CIBres/PBMC",x,".csv", sep=""),h=T)
# } ))
# colnames(PBMC.cib.ss) <- str_replace_all(colnames(PBMC.cib.ss),"\\.", " " )
# colnames(PBMC.cib.ss)[10] <- "T cells regulatory  (Tregs)" 
# rownames(PBMC.cib.ss) <- paste("PBMC",1:20,sep=" ")


##Process by MIXTURE
FL.robust <- MIXTURE(expressionMatrix = FL, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = 5L)

PBMC.robust <- MIXTURE(expressionMatrix = PBMC, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = 5L)

#LS fit
FL.abbas <- MIXTURE(expressionMatrix = FL, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = 3L)

PBMC.abbas <- MIXTURE(expressionMatrix = PBMC, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = 5L)
# #DTANGLE
# dt_FL <- dtangle(Y=newman_fl$data$log[1:14,] , reference = newman_fl$data$log[-c(1:14),])
# 
# dt_PBMC <- dtangle(Y=newman_pbmc$data$log[1:20,] , reference = newman_pbmc$data$log[-c(1:20),])
##RLM
abis_FL <- MIXTURE(expressionMatrix = FL, signatureMatrix =  LM22, functionMixture =  rlm.abis, useCores = 1L)

abis_PBMC <- MIXTURE(expressionMatrix = PBMC, signatureMatrix =  LM22, functionMixture =  rlm.abis, useCores = 1L)

quanti_FL <- out.quanti <-deconvolute_quantiseq.default(mix.mat = FL, 
                                                          arrays = FALSE, 
                                                          signame = "/home/elmer/Dropbox/IDEAS/cibersort/MIXTURE/Data/LM22", 
                                                          tumor = FALSE, 
                                                          mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")

quanti_PBMC <- out.quanti <-deconvolute_quantiseq.default(mix.mat = PBMC, 
                                                          arrays = FALSE, 
                                                          signame = "/home/elmer/Dropbox/IDEAS/cibersort/MIXTURE/Data/LM22", 
                                                          tumor = FALSE, 
                                                          mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
##Since FL data contains only B , CD8 and CD$ ct, we merge and summarise such CEll types according to Newman et al Supplementary File 
#https://www.nature.com/articles/nmeth.3337#supplementary-information
#
fl.composite <- c("B", "R","CD8T","CD4T","R","R","R","R","R","R","R","R")
FL.mixture <- CompositeMixture(newman_fl$annotation$mixture[,c(3,4,2,1,5,6,7,8,9,10,11,12)],fl.composite)[1:14,]
FL.mixture <- FL.mixture[,order(colnames(FL.mixture))]

FL.cs.c <- CompositeMixture(FL.cib.site, LM22.to.FLtypes)[,order(colnames(FL.mixture))]
# FL.cib.c.ss <-CompositeMixture(data.matrix(FL.cib.ss[,2:23]), cell.types.names11)
# cbind(colnames(FL.cib.site),colnames(FL.cib.ss[,2:23]))
FL.r.c <- CompositeMixture(GetMixture(FL.robust,"prop"), LM22.to.FLtypes)[,order(colnames(FL.mixture))]
FL.a.c <- CompositeMixture(GetMixture(FL.abbas,"prop"), LM22.to.FLtypes)[,order(colnames(FL.mixture))]
# FL.dt.c <- CompositeMixture(dt_FL$estimates, LM22.to.FLtypes)[,order(colnames(FL.mixture))]
# FL.dt.c.ss <- CompositeMixture(dt_FLss, cell.types.names11)

FL.abis <- CompositeMixture(GetMixture(abis_FL,"prop"), LM22.to.FLtypes)[,order(colnames(FL.mixture))]
FL.quanti <- CompositeMixture(quanti_FL[,-c(1,24)], LM22.to.FLtypes)[,order(colnames(FL.mixture))]


df.FL <- data.frame(p = c(as.numeric(FL.cs.c), as.numeric(FL.a.c), as.numeric(FL.r.c),as.numeric(FL.abis), as.numeric(FL.quanti)) ,
                    model = factor( rep(c("CIBERSORT", "ABBAS","MIXTURE","ABIS","QUANTISEQ"), 
                                        each = length(as.numeric(FL.cs.c))), 
                                    levels = c("ABBAS","ABIS","CIBERSORT","MIXTURE","QUANTISEQ")),
                    CT = rep(colnames(FL.mixture),each= nrow(FL.cs.c),times = 5),
                    truth = rep(as.numeric(FL.mixture), times = 5))

df.FL$dif <- df.FL$p - df.FL$truth 


wilcox.test(dif~model, subset(df.FL, truth > 0 & model %in% c("CIBERSORT", "MIXTURE")), paired =TRUE, alternative = "less")

# Wilcoxon signed rank test with continuity correction
# 
# data:  dif by model
# V = 2408, p-value = 0.9998
# alternative hypothesis: true location shift is less than 0
wilcox.test(dif~model, subset(df.FL, truth > 0 & model %in% c("ABIS", "MIXTURE")), paired =TRUE, alternative = "less")
# Wilcoxon signed rank test
# 
# data:  dif by model
# V = 399, p-value = 0.2597
# alternative hypothesis: true location shift is less than 0
wilcox.test(dif~model, subset(df.FL, truth > 0 & model %in% c("QUANTISEQ", "MIXTURE")), paired =TRUE, alternative = "less")

wilcox.test(dif~model, subset(df.FL, truth == 0 & model %in% c("CIBERSORT", "MIXTURE")), paired =TRUE, alternative = "greater")
# Wilcoxon signed rank test with continuity correction
# 
# data:  dif by model
# V = 2408, p-value = 0.0002182
# alternative hypothesis: true location shift is greater than 0
wilcox.test(dif~model, subset(df.FL, truth == 0 & model %in% c("QUANTISEQ", "MIXTURE")), paired =TRUE, alternative = "greater")


df.FL.noceros <- subset(df.FL, truth > 0)
df.sum.ceros <- ddply(subset(df.FL, truth == 0), .(model), summarise, 
      mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T), max = max(dif, na.rm=T), min= min(dif, na.rm=T), q3 = quantile(dif, 0.95),
      cor = cor(truth,p))
rownames(df.sum.ceros) <- df.sum.ceros[,1]
round(df.sum.ceros[,-1],3)
#true coefficientes == 0
#            mean    sd   max   min    q3 cor
# ABBAS     0.181 0.060 0.297 0.096 0.286  NA
# ABIS      0.086 0.044 0.193 0.029 0.157  NA
# CIBERSORT 0.170 0.046 0.288 0.106 0.242  NA
# MIXTURE   0.138 0.049 0.259 0.080 0.220  NA
# QUANTISEQ 0.058 0.076 0.259 0.000 0.201  NA
# 
df.FL.summary <- ddply(df.FL, .(model), summarise, min= min(dif,na.rm=T),mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T), max = max(dif,na.rm = T),
                       cor = cor(truth,p))
rownames(df.FL.summary) <- df.FL.summary[,1]
round(df.FL.summary[,-1],3)
#              min mean    sd   max   cor
# ABBAS     -0.378    0 0.186 0.297 0.836
# ABIS      -0.345    0 0.132 0.193 0.933
# CIBERSORT -0.384    0 0.174 0.288 0.901
# MIXTURE   -0.359    0 0.158 0.259 0.895
# QUANTISEQ -0.402    0 0.148 0.460 0.927
# 
df.FL.nc.summary <- ddply(df.FL.noceros, .(model), summarise, min= min(dif,na.rm=T),mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T), max = max(dif,na.rm = T),
                       cor = cor(truth,p))
rownames(df.FL.nc.summary) <- df.FL.nc.summary[,1]
round(df.FL.nc.summary[,-1],3)
#true coefficientes > 0
#              min   mean    sd   max   cor
# ABBAS     -0.378 -0.061 0.174 0.201 0.866
# ABIS      -0.345 -0.029 0.139 0.172 0.925
# CIBERSORT -0.384 -0.057 0.164 0.155 0.921
# MIXTURE   -0.359 -0.046 0.155 0.205 0.902
# QUANTISEQ -0.402 -0.019 0.162 0.460 0.936
# 

p.FL <- ggplot(df.FL, aes(truth, dif)) + 
  geom_point(na.rm=TRUE, aes(colour = CT, shape = CT)) +
   geom_hline(data=df.FL.nc.summary,aes(yintercept=round(mean,3)), col = "red") +
   geom_hline(data=df.FL.nc.summary,aes(yintercept=round(mean+2*sd,3)), col = "red", linetype = "dashed") + 
   geom_hline(data=df.FL.nc.summary,aes(yintercept=round(mean-2*sd,3)), col = "red", linetype = "dashed") +
  geom_smooth(data = df.FL.noceros, span = 0.5) + 
  labs( x = "True Flow Cytometry Derived (FCD) proportions",y = "Error")  +
  theme_bw()+  ggtitle("B")+
  facet_wrap(~model, nrow = 1)


##Generate the figures#####
setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/NuevoFL_BlandAltmanLM22.eps")##we can manage better
 print(p.FL)
dev.off()
# 
# png("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/NuevoFL_BlandAltman.png")##we can manage better
# print(p.FL)
# dev.off()

fl.cor<-ggplot(df.FL, aes(truth, p)) + 
  geom_point(na.rm=TRUE, aes(colour = CT, shape = CT)) +
  # geom_smooth(data = df.FL.noceros, span = 0.5) + 
  labs( x = "True Flow Cytometry Derived (FCD) proportions",y = "Error")  +
  theme_bw()+  ggtitle("B")+
  geom_abline(intercept = 0, slope = 1)+
  stat_smooth(method="lm", se=T, colour = "red",linetype="dashed", size=0.7)+
  facet_wrap(~model, nrow = 1)

setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/NuevoFL_LinearLM22.eps")##we can manage better
 print(fl.cor)
dev.off()

fl.bxp <- ggplot(df.FL, aes(model, dif, colour = model)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), trim = FALSE) + 
  geom_hline(yintercept = 0, colour = "red", linetype = "dotted") + labs( x = "True Flow Cytometry Derived (FCD) proportions",y = "Error")  +
  theme_bw()+  ggtitle("LM22 - Error distribution")

setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/NuevoFL_ErrorsDistributionLM22.eps")##we can manage better
 print(fl.bxp)
dev.off()
  


##PBMC Newman et al
##acomodar los nombres de las columnas para que sean iguales a los de LM22
library(stringr)
colnames(newman_pbmc$annotation$mixture) <- str_replace_all(colnames(newman_pbmc$annotation$mixture),"\\."," ")
colnames(newman_pbmc$annotation$mixture)
colnames(newman_pbmc$annotation$mixture)[which(str_detect(colnames(newman_pbmc$annotation$mixture), "T cells regulatory"))] <- "T cells regulatory (Tregs)"
newman_pbmc$annotation$mixture <- newman_pbmc$annotation$mixture[colnames(LM22)]

colnames(PBMC.cib.site)  <- str_replace_all( colnames(PBMC.cib.site), "\\." , " ")
colnames(PBMC.cib.site)[ which(str_detect(colnames(PBMC.cib.site), "T cells regulatory"))] <- "T cells regulatory (Tregs)"

cn <- colnames(data.matrix(newman_pbmc$annotation$mixture[1:20,]))
cn[newman_pbmc$annotation$mixture[1,]==0] <- "R"
cbind(colnames(LM22),colnames(data.matrix(newman_pbmc$annotation$mixture[1:20,])),cn)

PBMC.data <- CompositeMixture(data.matrix(newman_pbmc$annotation$mixture[1:20,]), cn)
colnames(PBMC.data)


PBMC.cs.c <-CompositeMixture(PBMC.cib.site, cn)
PBMC.r.c <- CompositeMixture(GetMixture(PBMC.robust,"prop"),cn)
PBMC.a.c <- CompositeMixture(GetMixture(PBMC.abbas,"prop"),cn)
# PBMC.dt.c <- dt_PBMC$estimates[,colnames(LM22)]
PBMC.abis <- CompositeMixture(GetMixture(abis_PBMC,"prop"),cn)
PBMC.quanti.e <- CompositeMixture(data.matrix(quanti_PBMC[,-c(1,24)]),cn)

df.t <- data.frame(x=rep(as.numeric(PBMC.data),5), 
                   y = c(as.numeric(PBMC.cs.c),as.numeric(PBMC.r.c),as.numeric(PBMC.a.c),
                         as.numeric(PBMC.abis), as.numeric(PBMC.quanti.e)),
                   model = factor(rep(c("CIBERSORT","MIXTURE","ABBAS","ABIS","QUANTISEQ"), each=length(as.numeric(PBMC.cs.c))),
                                  levels = c("ABBAS","ABIS","CIBERSORT","MIXTURE","QUANTISEQ")),
                   CT = rep(colnames(PBMC.data),each = ncol(PBMC.data),times = 5))



df.t$dif <- df.t$y-df.t$x
df.t.summary <- ddply(subset(df.t, x>0), .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
library(ggplot2)
df.PBMC <- df.t
df.PBMC.summary <- ddply(df.PBMC, .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
df.PBMC.summary <- ddply(subset(df.PBMC, x>0), .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
#true coefficientes > 0
#       model         mean         sd
# 1     ABBAS -0.047056061 0.12019068
# 2      ABIS -0.019395463 0.15119402
# 3   DTANGLE -0.029325333 0.08241455
# 4 CIBERSORT -0.020569976 0.09303382
# 5   MIXTURE -0.009712779 0.09799701
# 6   QUANTISEQ -0.032148524 0.16762496
ddply(subset(df.PBMC, x>0), .(model), summarise, cor = cor(x,y))
# model       cor
# 1     ABBAS 0.2739105
# 2      ABIS 0.3576994
# 3   DTANGLE 0.5641325
# 4 CIBERSORT 0.5733314
# 5   MIXTURE 0.6235775
# 6 QUANTISEQ 0.1751661
# 
ddply(subset(df.PBMC, x == 0), .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
#true coefficientes == 0
#       model        mean         sd
# 1     ABBAS 0.032577273 0.06197797
# 2      ABIS 0.013427628 0.06609093
# 3   DTANGLE 0.020302154 0.02710329
# 4 CIBERSORT 0.014240776 0.02894889
# 5   MIXTURE 0.006724231 0.02589453
# 6 QUANTISEQ 0.022256671 0.07770183

wilcox.test(dif~model, subset(df.PBMC, x > 0 & model %in% c("CIBERSORT", "MIXTURE")), paired =TRUE, alternative = "less")
# Wilcoxon signed rank test with continuity correction
# 
# data:  dif by model
# V = 4300, p-value = 0.007364
# alternative hypothesis: true location shift is less than 0
wilcox.test(dif~model, subset(df.PBMC, x > 0 & model %in% c("ABIS", "MIXTURE")), paired =TRUE, alternative = "less")
# Wilcoxon signed rank test with continuity correction
# 
# data:  dif by model
# V = 6529, p-value = 0.01051
# alternative hypothesis: true location shift is less than 0
wilcox.test(dif~model, subset(df.PBMC, x == 0 & model %in% c("CIBERSORT", "MIXTURE")), paired =TRUE, alternative = "greater")
# Wilcoxon signed rank test with continuity correction
# 
# data:  dif by model
# V = 5064, p-value = 4.368e-09
# alternative hypothesis: true location shift is greater than 0

col.cel.types <- c("chocolate1", "brown4", "black",
                   "tan1","green4", "green2", "lawngreen", "olivedrab3", "olivedrab", "chartreuse4",
                   "goldenrod","gold4","yellow","violetred","orangered1","red",
                   "plum4","plum","navy","mediumblue","cyan",
                   "grey28")
colores <- data.frame(CT=as.character(unique(df.PBMC$CT)), Colores=col.cel.types)
rownames(colores) <- colores[,1]

p.PBMC <-ggplot(df.PBMC, aes(x, dif)) + 
  geom_point(na.rm=TRUE, aes(colour = CT)) + 
  # scale_fill_manual(values  = as.character(colores[levels(df.t$CT),2]))+
  # theme(legend.position = "none") +
  geom_hline(data=df.PBMC.summary,aes(yintercept=round(mean,3)),col="red") +
  geom_hline(data=df.PBMC.summary,aes(yintercept=round(mean+2*sd,3)), col = "red", linetype = "dashed") + 
  geom_hline(data=df.PBMC.summary,aes(yintercept=round(mean-2*sd,3)), col = "red", linetype = "dashed") + 
  geom_smooth(data = subset(df.t, df.t$x > 0), span =.7)+
  labs( x = "True Flow Cytometry Derived (FCD) proportions",y = "Error") +
  theme_bw()+  ggtitle("PBMC - LM22")+
  facet_wrap(~model, nrow = 1)

setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/PBMC_BlandAltmanLM22.eps")##we can manage better
 print(p.PBMC)
dev.off()
# 
# png("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/NuevoFL_BlandAltman.png")##we can manage better
# print(p.FL)
# dev.off()

PBMC.cor<-ggplot(df.PBMC, aes(x, y)) + 
  geom_point(na.rm=TRUE, aes(colour = CT, shape = CT)) +
  # geom_smooth(data = df.FL.noceros, span = 0.5) + 
  labs( x = "True Flow Cytometry Derived (FCD) proportions",y = "Predicted proportions")  +
  theme_bw()+  ggtitle("PBMC - LM22 - Linear Regression")+
  geom_abline(intercept = 0, slope = 1)+
  stat_smooth(method="lm", se=T, colour = "red",linetype="dashed", size=0.7)+
  facet_wrap(~model, nrow = 1)

setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/NuevoPBMC_LinearLM22.eps")##we can manage better
 print(PBMC.cor)
dev.off()


## correlation
ddply(df.PBMC, .(model), summarise, cor = cor(x,y))
# model       cor
# 1     ABBAS 0.2875958
# 2      ABIS 0.4416744
# 3   DTANGLE 0.6688340
# 4 CIBERSORT 0.6720866
# 5   MIXTURE 0.7232487
# 6 QUANTISEQ 0.2647742
ddply(subset(df.PBMC,x>0) , .(model), summarise, cor = cor(x,y))


ddply(df.FL, .(model), summarise, cor = cor(truth,p))
# model       cor
# 1     ABBAS 0.8674040
# 2      ABIS 0.8913935
# 3   DTANGLE 0.8769220
# 4 CIBERSORT 0.8980177
# 5   MIXTURE 0.9023532
# 6    QUANTI 0.9514655
ddply(subset(df.FL, truth >0), .(model), summarise, cor = cor(truth,p))
# model       cor
# 1     ABBAS 0.9561770
# 2      ABIS 0.9456332
# 3   DTANGLE 0.9035471
# 4 CIBERSORT 0.9403999
# 5   MIXTURE 0.9425537
# 6    QUANTI 0.9665101

# png("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/NewPBMC_BlandAltman.png")##we can manage better
# print(p.PBMC)
# dev.off()
# 
# setEPS()
# postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/NewPBMC_BlandAltman.eps")##we can manage better
# print(p.PBMC)
# dev.off()






