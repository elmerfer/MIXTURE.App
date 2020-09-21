##REAL DATA from Noewman
rm(list=ls())
library(dtangle)
library(immunedeconv)
source('Utils/MIXTURE.DEBUG_V0.1.R')
source('Utils/dm_utils.R')

ROCAnalisis <- function(df.ba){
  res <- do.call(rbind,lapply(unique(df.ba$Method), function(x){
    test <- subset(df.ba, Method == x)
    ROCeval(True = factor(test$betasim>0, levels=c(TRUE,FALSE)), 
            Predicted  = factor(test$betahat>0, levels=c(TRUE,FALSE)))$Stats
  }))
  rownames(res) <- unique(df.ba$Method)
  return(data.frame(res))
}

setbetacol <- function(df){
  df$beta <- ifelse(df$betasim >0, ">0","=0")
  return(df[,order(colnames(df))])
}

BA.plot <- function(df, title = NULL, xlab = FALSE, plot=FALSE, scale = NULL){
  df.summary <- ddply(df, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
                      min= min(betahat,na.rm=T), max = max(betahat,na.rm=T), cor = cor(betahat, betasim))
  qq <-c(0.1,NA, 0.8*max(df$betasim))
  qq[2] <- qq[1]+ (qq[3]-qq[1])/2
  qq <- round(qq,2)
  p <- ggplot(df, aes(betasim, difs)) + 
    geom_point(size =0.7,aes(colour = CT)) + 
    geom_hline(yintercept = 0, colour = "blue",linetype = "dotted") +
    geom_hline(data=df.summary,aes(yintercept=round(mean,3)), color = "red", linetype = "solid") +
    geom_hline(data=df.summary,aes(yintercept=round(mean+2*sd,3)), color = "red", linetype = "dashed") + 
    geom_hline(data=df.summary,aes(yintercept=round(mean-2*sd,3)), color = "red", linetype = "dashed") +
    geom_smooth(aes(betasim, difs), data=subset(df,betasim>0)) +  scale_x_continuous(
      breaks = qq)+ theme_bw() + theme(axis.text.x = element_text(size=9),axis.text.y = element_text(size=9), legend.position = "none")+
    ylim(-.45, 0.6)
  if(!is.null(scale)) p <- p + scale
  if(!is.null(title)) p <- p + ggtitle(title)
  if(xlab) {
    p <- p + labs( x = "Simulated coefficients.",y = "Error")
  }else p <- p + labs( x = "",y = "Error")
  p <-  p +    facet_wrap(~Method, nrow = 1)
  if(plot) print(p)
  return(p)
}

Cor.plot <- function(df, title = NULL, xlab = FALSE, plot=FALSE, scale=NULL){
  df.summary <- ddply(df, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
                      min= min(betahat,na.rm=T), max = max(betahat,na.rm=T), cor = cor(betahat, betasim))
  qq <-c(0.1,NA, 0.8*max(df$betasim))
  qq[2] <- qq[1]+ (qq[3]-qq[1])/2
  qq <- round(qq,2)
  
  p <- ggplot(df, aes(betasim, betahat)) + 
    geom_point(size =0.7, aes(colour = CT)) + 
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    geom_smooth(method = "lm") + scale_x_continuous(
      breaks = qq)+ theme_bw() + theme(axis.text.x = element_text(size=9),axis.text.y = element_text(size=9), legend.position = "none")
  
  if(!is.null(scale)) p <- p +  scale
  
  if(!is.null(title)) p <- p + ggtitle(title)
  if(xlab) {
    p <- p + labs( x = "Simulated coefficients.",y = "Est. coeffs")
  }else p <- p + labs( x = "",y = "Est. coeffs")
  p <-  p +    facet_wrap(~Method, nrow = 1)
  
  if(plot) print(p)
  return(p)
}

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

FL.raw <- newman_fl$annotation$mixture[1:14,c(3,4,2,1,5,6,7,8,9,10,11,12)]
View(FL.raw)
colnames(FL.raw)
LM22.original.22 <-  c("Bn"     ,"Bm"      ,"Pcs","CD8T","CD4T","CD4Tr","CD4Ta","Tfh","Tregs","gdt","Nkr","Nka","MM","M0","M1","M2","Dcr","Dca","Mcr","Mca","E","N")
LM22.to.FLtypes <-    c("B", "B", "P",  "CD8T","CD4T","CD4T", "CD4T", "CD4T","CD4T","TGD","Nk","Nk","Mo","M","M","M", "D",  "D",  "Ma",  "Ma", "E","N")

cbind(LM22.original.22,LM22.to.FLtypes)

FL.ct.raw <- c("B","D","CD8T","CD4T","E","M","Ma","Mo","N","NK","P","TGD")
FL.composite <- CompositeMixture(newman_fl$annotation$mixture[,c(3,4,2,1,5,6,7,8,9,10,11,12)],FL.ct.raw)[1:14,]
colnames(FL.composite)

FL.composite <- FL.composite[,order(colnames(FL.composite))]

FL.cibersort <- CompositeMixture(FL.cib.site, LM22.to.FLtypes)
FL.cibersort <- FL.cibersort[, order(colnames(FL.cibersort))]
colnames(FL.cibersort)
# FL.cib.c.ss <-CompositeMixture(data.matrix(FL.cib.ss[,2:23]), cell.types.names11)
# cbind(colnames(FL.cib.site),colnames(FL.cib.ss[,2:23]))
FL.mixture <- CompositeMixture(GetMixture(FL.robust,"prop"), LM22.to.FLtypes)
FL.mixture <- FL.mixture[,order(colnames(FL.mixture))]
colnames(FL.mixture)
FL.ABB <- CompositeMixture(GetMixture(FL.abbas,"prop"), LM22.to.FLtypes)[,order(colnames(FL.composite))]

# FL.dt.c <- CompositeMixture(dt_FL$estimates, LM22.to.FLtypes)[,order(colnames(FL.mixture))]
# FL.dt.c.ss <- CompositeMixture(dt_FLss, cell.types.names11)

FL.ABI <- CompositeMixture(GetMixture(abis_FL,"prop"), LM22.to.FLtypes)[,order(colnames(FL.mixture))]
FL.QUA <- CompositeMixture(quanti_FL[,-c(1,24)], LM22.to.FLtypes)[,order(colnames(FL.mixture))]

LM22.FL <- data.frame(betahat = c(c(FL.cibersort), c(FL.ABB), c(FL.mixture),
                                c(FL.ABI), c(FL.QUA)) ,
                    Method = factor( rep(c("CIBERSORT", "ABBAS","MIXTURE","ABIS","QUANTISEQ"), 
                                         each = length(as.numeric(FL.composite))), 
                                     levels = c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")),
                    CT = rep(colnames(FL.composite),each= nrow(FL.composite),times = 5),
                    betasim = rep(as.numeric(FL.composite), times = 5))

LM22.FL$difs <- LM22.FL$betahat - LM22.FL$betasim
LM22.FL <- setbetacol(LM22.FL)
ROCAnalisis(LM22.FL)
#           TP  TN FP FN     Se    Sp   PPV    NPV    DO DOSeSp    ERG       F1
# CIBERSORT 42  47 79  0 100.00 37.30 34.71 100.00 90.52  62.70 0.6270 51.53292
# ABBAS     24  80 46 18  57.14 63.49 34.29  81.63 88.46  56.30 0.7937 42.85969
# MIXTURE   42  92 34  0 100.00 73.02 55.26 100.00 52.25  26.98 0.2698 71.18382
# ABIS      35  76 50  7  83.33 60.32 41.18  91.57 73.37  43.04 0.5635 55.12054
# QUANTISEQ 17 110 16 25  40.48 87.30 51.52  81.48 79.98  60.86 0.7222 45.33760





LM22.BA.FL <- BA.plot(LM22.FL,plot=T)
LM22.Cor.FL <- Cor.plot(LM22.FL,plot=T,scale = scale_y_continuous(
  breaks = c(-2,-1,0,1,2),
  label = c("-2.0","-1.0", "0.0", "1.0","2.0")))

ggplot(LM22.FL,aes(x=difs, group=beta,fill=beta)) + geom_density(alpha=0.4,size = 0.2) + labs(x="error", fill = expression(beta))+facet_grid(~ Method) 



df.FL.noceros <- subset(LM22.FL, betasim > 0)
df.sum.ceros <- ddply(subset(LM22.FL, betasim == 0), .(Method), summarise, 
      mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), max = max(difs, na.rm=T), min= min(difs, na.rm=T), q3 = quantile(difs, 0.95),
      cor = cor(betasim,betahat))
rownames(df.sum.ceros) <- df.sum.ceros[,1]
round(df.sum.ceros[,-1],3)
#true coefficientes == 0
#             mean    sd   max    min    q3 cor
# ABBAS     0.042 0.076 0.297  0.000 0.211  NA
# ABIS      0.032 0.072 0.241 -0.061 0.178  NA
# QUANTISEQ 0.011 0.039 0.259  0.000 0.090  NA
# CIBERSORT 0.025 0.038 0.205  0.000 0.098  NA
# MIXTURE   0.018 0.042 0.250  0.000 0.104  NA
# 
df.FL.summary <- ddply(LM22.FL, .(Method), summarise, min= min(difs,na.rm=T),mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), max = max(difs,na.rm = T),
                       cor = cor(betasim,betahat))
rownames(df.FL.summary) <- df.FL.summary[,1]
round(df.FL.summary[,-1],3)
#             min mean    sd   max   cor
# ABBAS     -0.427    0 0.115 0.297 0.861
# ABIS      -0.403    0 0.106 0.241 0.879
# QUANTISEQ -0.333    0 0.083 0.384 0.943
# CIBERSORT -0.442    0 0.090 0.205 0.954
# MIXTURE   -0.358    0 0.078 0.250 0.953
# 
df.FL.nc.summary <- ddply(subset(LM22.FL, betasim>0), .(Method), summarise, min= min(difs,na.rm=T),mean = mean(difs, na.rm=T), 
                          sd = sd(difs, na.rm=T), max = max(difs,na.rm = T),
                       cor = cor(betasim,betahat))
rownames(df.FL.nc.summary) <- df.FL.nc.summary[,1]
round(df.FL.nc.summary[,-1],3)
#true coefficientes > 0
#             min   mean    sd   max   cor
# ABBAS     -0.427 -0.125 0.121 0.103 0.941
# ABIS      -0.403 -0.096 0.134 0.131 0.920
# QUANTISEQ -0.333 -0.034 0.147 0.384 0.949
# CIBERSORT -0.442 -0.074 0.144 0.090 0.961
# MIXTURE   -0.358 -0.053 0.124 0.114 0.957
# 



##PBMC Newman et al
##acomodar los nombres de las columnas para que sean iguales a los de LM22
library(stringr)
colnames(newman_pbmc$annotation$mixture) <- str_replace_all(colnames(newman_pbmc$annotation$mixture),"\\."," ")
colnames(newman_pbmc$annotation$mixture)
colnames(newman_pbmc$annotation$mixture)[which(str_detect(colnames(newman_pbmc$annotation$mixture), "T cells regulatory"))] <- "T cells regulatory (Tregs)"
newman_pbmc$annotation$mixture <- newman_pbmc$annotation$mixture[colnames(LM22)]

colnames(PBMC.cib.site)  <- str_replace_all( colnames(PBMC.cib.site), "\\." , " ")
colnames(PBMC.cib.site)[ which(str_detect(colnames(PBMC.cib.site), "T cells regulatory"))] <- "T cells regulatory (Tregs)"


PBMC.composite <- data.matrix(newman_pbmc$annotation$mixture[1:20,])
colnames(PBMC.composite)
View(PBMC.composite)


PBMC.cibersort <-PBMC.cib.site
PBMC.mixture <- GetMixture(PBMC.robust,"prop")
PBMC.ABB <- GetMixture(PBMC.abbas,"prop")
# PBMC.dt.c <- dt_PBMC$estimates[,colnames(LM22)]
PBMC.ABI <- GetMixture(abis_PBMC,"prop")
PBMC.QUA <- data.matrix(quanti_PBMC[,-c(1,24)])
colnames(PBMC.QUA) <- str_replace_all(colnames(PBMC.QUA),"\\."," ")

cbind(colnames(PBMC.composite),colnames(PBMC.cibersort),colnames(PBMC.mixture),colnames(PBMC.ABB),colnames(PBMC.ABI),colnames(PBMC.QUA))

LM22.PBMC <- data.frame(betahat = c(c(PBMC.cibersort), c(PBMC.ABB), c(PBMC.mixture),
                               c(PBMC.ABI), c(PBMC.QUA)) ,
                   Method = factor( rep(c("CIBERSORT", "ABBAS","MIXTURE","ABIS","QUANTISEQ"), 
                                        each = length(as.numeric(PBMC.composite))), 
                                    levels = c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")),
                   CT = rep(colnames(PBMC.composite),each= nrow(PBMC.composite),times = 5),
                   betasim = rep(as.numeric(PBMC.composite), times = 5))

LM22.PBMC$difs <- LM22.PBMC$betahat - LM22.PBMC$betasim
LM22.PBMC <- setbetacol(LM22.PBMC)

ROCAnalisis(LM22.PBMC)
#           TP  TN  FP  FN    Se    Sp   PPV   NPV    DO DOSeSp    ERG       F1
# CIBERSORT 148 149 111  32 82.22 57.31 57.14 82.32 65.48  46.24 0.6047 67.42323
# ABBAS      89 164  96  91 49.44 63.08 48.11 64.31 88.80  62.61 0.8748 48.76593
# MIXTURE   129 223  37  51 71.67 85.77 77.71 81.39 42.99  31.70 0.4256 74.56789
# ABIS      129 141 119  51 71.67 54.23 52.02 73.44 76.84  53.83 0.7410 60.28415
# QUANTISEQ  73 225  35 107 40.56 86.54 67.59 67.77 76.18  60.94 0.7290 50.69719


LM22.BA.PBMC <- BA.plot(LM22.PBMC,plot=T)
LM22.Cor.PBMC <- Cor.plot(LM22.PBMC,plot=T,scale = scale_y_continuous(
  breaks = c(-2,-1,0,1,2),
  label = c("-2.0","-1.0", "0.0", "1.0","2.0")))

ggplot(LM22.PBMC,aes(x=difs, group=beta,fill=beta)) + geom_density(alpha=0.4,size = 0.2) + labs(x="error", fill = expression(beta))+facet_grid(~ Method) 



df.t$dif <- df.t$y-df.t$x
df.t.summary <- ddply(subset(df.t, x>0), .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
library(ggplot2)

df.PBMC <- df.t
colnames(df.t) <-  c("betasim","betahat","Method","CT","dif")
ROCAnalisis(df.t)
#           TP TN FP  FN    Se Sp   PPV  NPV     DO DOSeSp    ERG       F1
# CIBERSORT 148  0 20  32 82.22  0 88.10 0.00 143.03 101.57 1.1778 85.05850
# MIXTURE   129  3 17  51 71.67 15 88.36 5.56 130.70  89.60 1.1333 79.14468
# ABBAS      89  0 20  91 49.44  0 81.65 0.00 151.30 112.05 1.5056 61.58786
# ABIS      129  0 20  51 71.67  0 86.58 0.00 144.85 103.94 1.2833 78.42260
# QUANTISEQ  73  0 20 107 40.56  0 78.49 0.00 154.91 116.33 1.5944 53.48264
df.PBMC.summary <- ddply(df.PBMC, .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
df.PBMC.summary <- ddply(subset(df.PBMC, x>0), .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
df.PBMC.summary
#true coefficientes > 0
# model         mean         sd
# 1     ABBAS -0.047056061 0.12019068
# 2      ABIS -0.019395463 0.15119402
# 3 CIBERSORT -0.020569976 0.09303382
# 4   MIXTURE -0.009712779 0.09799701
# 5 QUANTISEQ -0.032148524 0.16762496
ddply(subset(df.PBMC, x>0), .(model), summarise, cor = cor(x,y))
# model       cor
# 1     ABBAS 0.2739105
# 2      ABIS 0.3576994
# 3 CIBERSORT 0.5733314
# 4   MIXTURE 0.6235775
# 5 QUANTISEQ 0.1751661
# 
ddply(subset(df.PBMC, x == 0), .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
#true coefficientes == 0
# model       mean         sd
# 1     ABBAS 0.42350455 0.06928903
# 2      ABIS 0.17455916 0.07456590
# 3 CIBERSORT 0.18513008 0.06023039
# 4   MIXTURE 0.08741501 0.08035369
# 5 QUANTISEQ 0.28933672 0.16170951

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
colnames(df.PBMC)



## correlation
ddply(df.PBMC, .(model), summarise, cor = cor(x,y))
# model         cor
# 1     ABBAS -0.08104465
# 2      ABIS  0.27166788
# 3 CIBERSORT  0.41641651
# 4   MIXTURE  0.58497915
# 5 QUANTISEQ  0.01545036

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
# 1     ABBAS 0.8356833
# 2      ABIS 0.9326970
# 3 CIBERSORT 0.9011704
# 4   MIXTURE 0.8954254
# 5 QUANTISEQ 0.9271279

save(LM22.FL,LM22.PBMC, file = "flow_cyotmetry_output_LM22.RData"  )





