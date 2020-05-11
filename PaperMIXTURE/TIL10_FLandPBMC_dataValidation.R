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
#"B" "CD8" "CD4" "NK"  "Mo"  "M0"  "M1"  "M2"  "D"   "M"   "E"   "N" 
#Representacion de los tipos celulares para TIL10
# colnames(TIL10)
# [1] "B.cells"         "Macrophages.M1"  "Macrophages.M2"  "Monocytes"       "Neutrophils"     "NK.cells"        "T.cells.CD4"    
# [8] "T.cells.CD8"     "Tregs"           "Dendritic.cells"

TIL10.ct <-   c("B","M1","M2","Mo","N","NK","CD4T","CD8T","Tregs","D")
##colnames for mapping TIL10 data to FL data
TIL10.ct.2 <- c("B","R" ,"R" ,"R" ,"R" ,"R","CD4T","CD8T","CD4T","R")
# 1] "T.CD4"         "T.CD8"         "B"             "Dendritic"     "Eosinophils"   "Macrophages"   "Mast"          "Monocytes"
# [9] "Neutrophils"   "NK"            "Plasma"        "T.Gamma.Delta"
fl.composite <- c("CD4T", "CD8T","B","R","R","R","R","R","R","R","R","R")

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
load("Data/TIL10.RData")

path.2.TIL10 <- system.file("extdata", "quantiseq", paste0("TIL10", "_signature.txt"), package = "immunedeconv", mustWork = TRUE)
#load the signature TIL10
TIL10 <- read.table(path.2.TIL10, header = TRUE, sep = "\t", row.names = 1)
rownames(TIL10)

save(TIL10, "./Data/TIL10.RData")

##prepare data for CIBERSORT and MIXTURE
FL <- t(2^newman_fl$data$log)[,1:14]
PBMC <- t(2^newman_pbmc$data$log)[,1:20]

FL_beta.table <- newman_fl$annotation$mixture[1:14,]
PBMC_beta.table <- newman_pbmc$annotation$mixture[1:20,]
colnames(FL_beta.table)[1:3] <- c("T.cells.CD4","T.cells.CD8","B.cells")
colnames(PBMC_beta.table)
colnames(TIL10)
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
FL.cib.site <- ReadCibersortWebResults(file = "./Data/CIBERSORT.Output_FL.TIL10.csv", type = "csv", nct =10)
PBMC.cib.site <- ReadCibersortWebResults(file = "Data/CIBERSORT.Output_PBMC.TIL10.csv", type = "csv", nct =10)

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
FL.robust <- MIXTURE(expressionMatrix = FL, signatureMatrix =  TIL10, functionMixture =  nu.svm.robust.RFE, useCores = 5L)

PBMC.robust <- MIXTURE(expressionMatrix = PBMC, signatureMatrix =  TIL10, functionMixture =  nu.svm.robust.RFE, useCores = 5L)

#LS fit
FL.abbas <- MIXTURE(expressionMatrix = FL, signatureMatrix =  TIL10, functionMixture =  ls.rfe.abbas, useCores = 3L)

PBMC.abbas <- MIXTURE(expressionMatrix = PBMC, signatureMatrix =  TIL10, functionMixture =  ls.rfe.abbas, useCores = 5L)
#DTANGLE
# dt_FL <- dtangle(Y=newman_fl$data$log[1:14,] , reference = newman_fl$data$log[-c(1:14),])
# 
# dt_PBMC <- dtangle(Y=newman_pbmc$data$log[1:20,] , reference = newman_pbmc$data$log[-c(1:20),])
##RLM
abis_FL <- MIXTURE(expressionMatrix = FL, signatureMatrix =  TIL10, functionMixture =  rlm.abis, useCores = 1L)

abis_PBMC <- MIXTURE(expressionMatrix = PBMC, signatureMatrix =  TIL10, functionMixture =  rlm.abis, useCores = 1L)

quanti_FL <- deconvolute_quantiseq.default(mix.mat = FL, 
                                                          arrays = FALSE, 
                                                          signame = "TIL10", 
                                                          tumor = FALSE, 
                                                          mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")

quanti_PBMC  <-deconvolute_quantiseq.default(mix.mat = PBMC, 
                                                          arrays = FALSE, 
                                                          signame = "TIL10", 
                                                          tumor = FALSE, 
                                                          mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")


##Since FL data contains only B , CD8 and CD$ ct, we merge and summarise such CEll types according to Newman et al Supplementary File 
#https://www.nature.com/articles/nmeth.3337#supplementary-information
FL.cs.c <- CompositeMixture(FL.cib.site, TIL10.ct.2 )
# FL.cib.c.ss <-CompositeMixture(data.matrix(FL.cib.ss[,2:23]), cell.types.names11)
# cbind(colnames(FL.cib.site),colnames(FL.cib.ss[,2:23]))
FL.r.c <- CompositeMixture(GetMixture(FL.robust,"prop"), TIL10.ct.2 )
FL.a.c <- CompositeMixture(GetMixture(FL.abbas,"prop"), TIL10.ct.2 )
# FL.dt.c <- CompositeMixture(dt_FL$estimates, TIL10.ct.2 )
# FL.dt.c.ss <- CompositeMixture(dt_FLss, cell.types.names11)

FL.abis <- CompositeMixture(GetMixture(abis_FL,"prop"), TIL10.ct.2 )
FL.quanti <- data.matrix(CompositeMixture(quanti_FL[,-c(1,11)], TIL10.ct.2 ))




colnames(newman_fl$annotation$mixture[1:14,c(3,4,2,1,5,6,7,8,9,10,11,12)])
##Mapping the cell types present in the FL data set 
fl.composite <- c("B", "R","CD8T","CD4T","R","R","R","R","R","R","R","R")
##The proportion matrix of FL data
FL.mixture <- CompositeMixture(newman_fl$annotation$mixture[,c(3,4,2,1,5,6,7,8,9,10,11,12)],fl.composite)[1:14,]

summary(FL.cs.c)
df.FL <- data.frame(betahat = c(as.numeric(FL.cs.c), as.numeric(FL.a.c), as.numeric(FL.r.c),
                                as.numeric(FL.abis), as.numeric(FL.quanti)) ,
                    Method = factor( rep(c("CIBERSORT", "ABBAS","MIXTURE","ABIS","QUANTISEQ"), 
                                        each = length(as.numeric(FL.cs.c))), 
                                    levels = c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")),
                    CT = rep(colnames(FL.mixture),each= nrow(FL.cs.c),times = 5),
                    betasim = rep(as.numeric(FL.mixture), times = 5))


df.FL$difs <- df.FL$betahat - df.FL$betasim 
df.FL$CT <- as.character(df.FL$CT)
df.FL$CT[df.FL$CT == "B"] <- "B cells"
df.FL$CT[df.FL$CT == "CD4T"] <- "T cells CD4"
df.FL$CT[df.FL$CT == "CD8T"] <- "T cells CD8"
df.FL$CT[df.FL$CT == "R"] <- "Other"
df.FL$CT <- factor(df.FL$CT, levels = c("B cells","T cells CD4","T cells CD8","Other"))

wilcox.test(difs~Method, subset(df.FL, betasim > 0 & Method %in% c("CIBERSORT", "MIXTURE")), paired =TRUE, alternative = "less")

# Wilcoxon signed rank test with continuity correction
# 
# data:  difs by Method
# V = 85, p-value = 6.453e-06
# alternative hypothesis: true location shift is less than 0
wilcox.test(difs~Method, subset(df.FL, betasim > 0 & Method %in% c("QUANTISEQ", "MIXTURE")), paired =TRUE, alternative = "less")

# Wilcoxon signed rank test
# 
# data:  dif by model
# V = 399, p-value = 0.2597
# alternative hypothesis: true location shift is less than 0
wilcox.test(dif~model, subset(df.FL, truth == 0 & model %in% c("CIBERSORT", "MIXTURE")), paired =TRUE, alternative = "greater")
# Wilcoxon signed rank test with continuity correction
# 
# data:  dif by model
# V = 2408, p-value = 0.0002182
# alternative hypothesis: true location shift is greater than 0

# 

df.FL.nc.summary <- ddply(subset(df.FL, betasim >0), .(Method), summarise, min = min(difs, na.rm=T), mean = mean(difs, na.rm=T), 
                          sd = sd(difs, na.rm=T), max = max(difs,na.rm=T), cor = cor(betasim,betahat))
rownames(df.FL.nc.summary) <- df.FL.nc.summary[,1]
round(df.FL.nc.summary[,-1],3)
#             min   mean    sd   max   cor
# ABBAS     -0.402 -0.033 0.211 0.328 0.785
# ABIS      -0.498  0.131 0.410 1.072 0.477
# QUANTISEQ -0.402  0.000 0.188 0.460 0.949
# CIBERSORT -0.384 -0.090 0.155 0.218 0.897
# MIXTURE   -0.352 -0.048 0.154 0.252 0.893

df.FL.nc.summary <- ddply(df.FL, .(Method), summarise, min = min(difs, na.rm=T), mean = mean(difs, na.rm=T), 
                          sd = sd(difs, na.rm=T), max = max(difs,na.rm=T), cor = cor(betasim,betahat))
rownames(df.FL.nc.summary) <- df.FL.nc.summary[,1]
round(df.FL.nc.summary[,-1],3)
# min mean    sd   max   cor
# ABBAS     -0.402    0 0.196 0.328 0.804
# ABIS      -1.528    0 0.458 1.072 0.584
# QUANTISEQ -0.402    0 0.162 0.460 0.950
# CIBERSORT -0.384    0 0.210 0.390 0.771
# MIXTURE   -0.352    0 0.164 0.312 0.866

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
    geom_smooth() +  scale_x_continuous(
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

BA.FL <- BA.plot(df.FL,plot=T)
Cor.FL <- Cor.plot(df.FL,plot=T,scale = scale_y_continuous(
  breaks = c(-2,-1,0,1,2),
  label = c("-2.0","-1.0", "0.0", "1.0","2.0")))



##PBMC Newman et al
##acomodar los nombres de las columnas para que sean iguales a los de LM22
library(stringr)
colnames(newman_pbmc$annotation$mixture) <- str_replace_all(colnames(newman_pbmc$annotation$mixture),"\\."," ")
colnames(newman_pbmc$annotation$mixture)
colnames(newman_pbmc$annotation$mixture)[which(str_detect(colnames(newman_pbmc$annotation$mixture), "T cells regulatory"))] <- "T cells regulatory (Tregs)"
newman_pbmc$annotation$mixture <- newman_pbmc$annotation$mixture[colnames(LM22)]

#Newman PBMC data only have the following FCD cells
colnames(newman_pbmc$annotation$mixture[1:20,])[newman_pbmc$annotation$mixture[1,]>0]
# [1] "B cells naive"                "B cells memory"               "T cells CD8"                  "T cells CD4 naive"            "T cells CD4 memory resting"  
# [6] "T cells CD4 memory activated" "T cells gamma delta"          "NK cells activated"           "Monocytes"  
#Thus we have to accommodate for that
lm22.to.til10 <- c("B cells","B cells","B cells","T cells CD8","T cells CD4","T cells CD4","T cells CD4",
                   "T cells CD4","Tregs","T cells CD8","NK cells","NK cells","Monocytes","Macrophages M1",
                   "Macrophages M1","Macrophages M2","Dendritic cells","Dendritic cells","R","R","R","Neutrophils") 
cbind(colnames(LM22),lm22.to.til10)
cn <- lm22.to.til10
cn[newman_pbmc$annotation$mixture[1,] == 0] <- "R"
cbind(colnames(newman_pbmc$annotation$mixture[1:20,]),cn)
pbmc.til10 <- unique(cn)

# [1,] "B cells"     "B cells naive"               
# [2,] "B cells"     "B cells memory"              
# [3,] "T cells CD8" "T cells CD8"                 
# [4,] "T cells CD4" "T cells CD4 naive"           
# [5,] "T cells CD4" "T cells CD4 memory resting"  
# [6,] "T cells CD4" "T cells CD4 memory activated"
# [7,] "T cells CD8" "T cells gamma delta"         
# [8,] "NK cells"    "NK cells activated"          
# [9,] "Monocytes"   "Monocytes"  




PBMC.quanti <- data.matrix(quanti_PBMC[,-c(1,24)])
colnames(quanti_PBMC)
# "B.cells"         "Macrophages.M1"  "Macrophages.M2"  "Monocytes"       "Neutrophils"     "NK.cells"        "T.cells.CD4"     "T.cells.CD8"     "Tregs"          
# [11] "Dendritic.cells"
# [1] "B cells"     "Monocytes"   "NK cells"    "R"           "T cells CD4" "T cells CD8"
#  "B cells"    "R"    "R"    "Monocytes"    "R"    "R"    "T cells CD4" "T cells CD8" "R" "R"
TIL10.ct.2 <- c("B cells","R","R","Monocytes","R","NK cells","T cells CD4","T cells CD8","R","R")
cbind(colnames(TIL10),TIL10.ct.2)

PBMC.data <- CompositeMixture(newman_pbmc$annotation$mixture[1:20,], cn)
PBMC.data <- PBMC.data[,order(colnames(PBMC.data))]
colnames(PBMC.data)
colnames(PBMC.cib.site)
PBMC.cs.c <-CompositeMixture(PBMC.cib.site,TIL10.ct.2)[,order(colnames(PBMC.data))]

PBMC.r.c <- CompositeMixture(GetMixture(PBMC.robust,"prop"),TIL10.ct.2)[,order(colnames(PBMC.data))]
PBMC.a.c <- CompositeMixture(GetMixture(PBMC.abbas,"prop"),TIL10.ct.2)[,order(colnames(PBMC.data))]
# PBMC.dt.c <- dt_PBMC$estimates[,colnames(LM22)]
PBMC.abis <- CompositeMixture(GetMixture(abis_PBMC,"prop"),TIL10.ct.2)[,order(colnames(PBMC.data))]
length(as.numeric(PBMC.abis))
PBMC.quanti.e <- CompositeMixture(data.matrix(quanti_PBMC[,-c(1,12)]), TIL10.ct.2)[,order(colnames(PBMC.data))]



df.t <- data.frame(betasim = rep(as.numeric(PBMC.data),5) , 
                   betahat = c(as.numeric(PBMC.cs.c),as.numeric(PBMC.r.c),as.numeric(PBMC.a.c),
                         as.numeric(PBMC.abis), as.numeric(PBMC.quanti.e)),
                   Method = factor(rep(c("CIBERSORT","MIXTURE","ABBAS","ABIS","QUANTISEQ"), each=length(as.numeric(PBMC.cs.c))),
                                  levels = c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")),
                   CT = rep(colnames(PBMC.data),each= nrow(PBMC.data), times = 5))

df.t$difs <- df.t$betahat-df.t$betasim
df.pbmc <- ddply(df.t, .(Method), summarise, min = min(difs, na.rm=T), mean = mean(difs, na.rm=T), 
                          sd = sd(difs, na.rm=T), max = max(difs,na.rm=T), cor = cor(betasim,betahat))
rownames(df.pbmc) <- df.pbmc[,1]
round(df.pbmc[,-1],3)
#               min mean    sd   max   cor
# ABBAS     -0.341    0 0.149 0.471 0.492
# ABIS      -0.624    0 0.202 0.722 0.466
# QUANTISEQ -0.620    0 0.354 0.955 0.039
# CIBERSORT -0.440    0 0.191 0.616 0.278
# MIXTURE1   -0.431    0 0.175 0.433 0.383

df.pbmc <- ddply(subset(df.t, betahat >0), .(Method), summarise, min = min(difs, na.rm=T), mean = mean(difs, na.rm=T), 
                 sd = sd(difs, na.rm=T), max = max(difs,na.rm=T), cor = cor(betasim,betahat))
rownames(df.pbmc) <- df.pbmc[,1]
round(df.pbmc[,-1],3)
#               min  mean    sd   max    cor
# ABBAS     -0.341 0.010 0.148 0.471  0.498
# ABIS      -0.366 0.020 0.194 0.722  0.495
# QUANTISEQ -0.221 0.439 0.367 0.955 -0.123
# CIBERSORT -0.440 0.005 0.197 0.616  0.231
# MIXTURE   -0.431 0.008 0.183 0.433  0.324


##Plot FIGURE.3
BA.FL <- BA.plot(df.FL,title = "A",plot=T)
Cor.FL <- Cor.plot(df.FL,title = "B",plot=T, scale=scale_y_continuous(
  breaks = c(-2,-1,0,1),
  label = c("-2.00","-1.00", "0.00", "1.00")))
BA.PBMC <- BA.plot(df.t,title = "C",plot=F)
Cor.PBMC <- Cor.plot(df.t,title = "D",plot=F, xlab = TRUE, scale=scale_y_continuous(
  breaks = c(0,0.5,1),
  label = c( "0.00", "0.50","1.00"))) 

pdf("FIG3.pdf",paper="a4", width = 10, height = 10)
grid.arrange(BA.FL, Cor.FL, BA.PBMC, Cor.PBMC, nrow = 4)
dev.off()
setEPS()

grid.arrange(BA.FL, Cor.FL, BA.PBMC, Cor.PBMC, nrow = 4)


postscript("FIG3.eps")
grid.arrange(BA.FL, Cor.FL, BA.PBMC, Cor.PBMC, nrow = 4)
dev.off()


rm(df.t.summary)
df.t.summary <- ddply(subset(df.t, betasim>0), .(Method), summarise, min = min(difs, na.rm=T), 
                       mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), max = max(difs,na.rm=T), cor = cor(betasim,betahat))
rownames(df.t.summary) <- df.t.summary[,1]
round(df.t.summary[,-1],3)
#             min   mean    sd   max    cor
# ABBAS     -0.341 -0.017 0.152 0.471  0.452
# ABIS      -0.624 -0.020 0.214 0.722  0.466
# QUANTISEQ -0.620  0.000 0.388 0.955 -0.086
# CIBERSORT -0.440 -0.018 0.200 0.616  0.214
# MIXTURE   -0.431 -0.020 0.179 0.404  0.353
##only FCD

df.t.summary <- ddply(df.t, betasim, .(Method), summarise, min = min(difs, na.rm=T), 
                      mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), max = max(difs,na.rm=T), cor = cor(betasim,betahat))
rownames(df.t.summary) <- df.t.summary[,1]
round(df.t.summary[,-1],3)
# min   mean    sd   max    cor
# ABBAS     -0.341 -0.017 0.152 0.471  0.452
# ABIS      -0.624 -0.020 0.214 0.722  0.466
# QUANTISEQ -0.620  0.000 0.388 0.955 -0.086
# CIBERSORT -0.440 -0.018 0.200 0.616  0.214
# MIXTURE   -0.431 -0.020 0.179 0.404  0.353

df.fl.til10 <- df.FL
df.pbmc.til10 <- df.t
save(df.fl.til10, df.pbmc.til10, file = "FL-PBMC_TIL10.RData")
