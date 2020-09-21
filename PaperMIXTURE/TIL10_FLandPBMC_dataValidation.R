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
 colnames(TIL10)
# [1] "B.cells"         "Macrophages.M1"  "Macrophages.M2"  "Monocytes"       "Neutrophils"     "NK.cells"        "T.cells.CD4"    
# [8] "T.cells.CD8"     "Tregs"           "Dendritic.cells"

TIL10.ct <-   c("B","M1","M2","Mo","N","NK","CD4T","CD8T","Tregs","D")
# TIL10.ct.inFL <-   c("B","M1","M2","Mo","N","NK","CD4T","CD8T","CD4T","D")
##colnames for mapping TIL10 data to FL data
#Para Folicular Lymphoma we only have B CD8 and CD4
TIL10.ct.2 <- c("B","R" ,"R" ,"R" ,"R" ,"R","CD4T","CD8T","CD4T","R")
summary(newman_fl$annotation$mixture[1:14,c(3,4,2,1,5,6,7,8,9,10,11,12)])
# B            Dendritic     T.CD8             T.CD4          Eosinophils  Macrophages      Mast     Monocytes  Neutrophils       NK        Plasma  T.Gamma.Delta
# Min.   :0.5402   Min.   :0   Min.   :0.01020   Min.   :0.05050   Min.   :0    Min.   :0    Min.   :0   Min.   :0   Min.   :0    Min.   :0   Min.   :0   Min.   :0    
# 1st Qu.:0.7219   1st Qu.:0   1st Qu.:0.02025   1st Qu.:0.09488   1st Qu.:0    1st Qu.:0    1st Qu.:0   1st Qu.:0   1st Qu.:0    1st Qu.:0   1st Qu.:0   1st Qu.:0    
# Median :0.8188   Median :0   Median :0.02940   Median :0.15515   Median :0    Median :0    Median :0   Median :0   Median :0    Median :0   Median :0   Median :0    
# Mean   :0.7852   Mean   :0   Mean   :0.03648   Mean   :0.17830   Mean   :0    Mean   :0    Mean   :0   Mean   :0   Mean   :0    Mean   :0   Mean   :0   Mean   :0    
# 3rd Qu.:0.8803   3rd Qu.:0   3rd Qu.:0.05750   3rd Qu.:0.22210   3rd Qu.:0    3rd Qu.:0    3rd Qu.:0   3rd Qu.:0   3rd Qu.:0    3rd Qu.:0   3rd Qu.:0   3rd Qu.:0    
# Max.   :0.9293   Max.   :0   Max.   :0.07140   Max.   :0.40230   Max.   :0    Max.   :0    Max.   :0   Max.   :0   Max.   :0    Max.   :0   Max.   :0   Max.   :0  
fl.composite

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
# TIL10 <- read.table(path.2.TIL10, header = TRUE, sep = "\t", row.names = 1)
# rownames(TIL10)

# save(TIL10, "./Data/TIL10.RData")
load("./Data/TIL10.RData")
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


##Process by MIXTURE TIL10
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


FL.raw <- newman_fl$annotation$mixture[1:14,c(3,4,2,1,5,6,7,8,9,10,11,12)]
View(FL.raw)
colnames(FL.raw)
# B"             "Dendritic"     "T.CD8"         "T.CD4"         "Eosinophils"   "Macrophages"   "Mast"          "Monocytes"     "Neutrophils"   "NK"           
# [11] "Plasma"        "T.Gamma.Delta"
FL.ct.raw <- c("B","D","CD8T","CD4T","E","M","Ma","Mo","N","NK","P","TGD")
colnames(TIL10)
# [1] "B.cells"         "Macrophages.M1"  "Macrophages.M2"  "Monocytes"       "Neutrophils"     "NK.cells"        "T.cells.CD4"     "T.cells.CD8"     "Tregs"          
# [10] "Dendritic.cells"
TIL10.ct.inFL <-   c("B","M","M","Mo","N","NK","CD4T","CD8T","CD4T","D")

sort(FL.ct.raw)
sort(TIL10.ct.inFL)


FL.composite <- CompositeMixture(FL.raw,FL.ct.raw)[,unique(TIL10.ct.inFL)]
FL.cibersort <- CompositeMixture(FL.cib.site,TIL10.ct.inFL)
FL.mixture <- CompositeMixture(GetMixture(FL.robust),TIL10.ct.inFL)
FL.ABB <- CompositeMixture(GetMixture(FL.abbas),TIL10.ct.inFL)
FL.ABI <- CompositeMixture(GetMixture(abis_FL),TIL10.ct.inFL)
FL.QUA <- CompositeMixture(quanti_FL[,-c(1,12)],TIL10.ct.inFL)




TIL10.FL <- data.frame(betahat = c(c(FL.cibersort), c(FL.ABB), c(FL.mixture),
                                   c(FL.ABI), c(FL.QUA)) ,
                       Method = factor( rep(c("CIBERSORT", "ABBAS","MIXTURE","ABIS","QUANTISEQ"), 
                                            each = length(as.numeric(FL.composite))), 
                                        levels = c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")),
                       CT = rep(colnames(FL.composite),each= nrow(FL.composite),times = 5),
                       betasim = rep(as.numeric(FL.composite), times = 5))

TIL10.FL$difs <- TIL10.FL$betahat - TIL10.FL$betasim
TIL10.FL <- setbetacol(TIL10.FL)
ROCAnalisis(TIL10.FL)
#           TP TN FP FN    Se     Sp    PPV   NPV    DO DOSeSp    ERG       F1
# CIBERSORT 39 38 32  3 92.86  54.29  54.93 92.68 65.00  46.26 0.5285 69.02767
# ABBAS     29 57 13 13 69.05  81.43  69.05 81.43 51.04  36.09 0.4952 69.05000
# MIXTURE   37 57 13  5 88.10  81.43  74.00 91.94 35.03  22.06 0.3047 80.43677
# ABIS      30 54 16 12 71.43  77.14  65.22 81.82 53.66  36.59 0.5143 68.18389
# QUANTISEQ 14 70  0 28 33.33 100.00 100.00 71.43 72.53  66.67 0.6667 49.99625

BA.FL <- BA.plot(TIL10.FL,plot=T)
Cor.FL <- Cor.plot(df.FL,plot=T,scale = scale_y_continuous(
  breaks = c(-2,-1,0,1,2),
  label = c("-2.0","-1.0", "0.0", "1.0","2.0")))

ggplot(TIL10.FL,aes(x=difs, group=beta,fill=beta)) + geom_density(alpha=0.4,size = 0.2) + labs(x="error", fill = expression(beta))+facet_grid(~ Method) 


dim(FL.cib.site)
View(newman_fl$annotation$mixture)

## PBMC
library(stringr)
colnames(newman_pbmc$annotation$mixture) <- str_replace_all(colnames(newman_pbmc$annotation$mixture),"\\."," ")
colnames(newman_pbmc$annotation$mixture)
# [1] "B cells naive"                "B cells memory"               "T cells CD8"                  "T cells CD4 naive"            "T cells CD4 memory resting"  
# [6] "T cells CD4 memory activated" "T cells gamma delta"          "NK cells activated"           "Monocytes"                    "Dendritic cells activated"   
# [11] "Dendritic cells resting"      "Eosinophils"                  "Macrophages M0"               "Macrophages M1"               "Macrophages M2"              
# [16] "Mast cells activated"         "Mast cells resting"           "Neutrophils"                  "NK cells resting"             "Plasma cells"                
# [21] "T cells follicular helper"    "T cells regulatory  Tregs " 
sum(colnames(LM22) %in% colnames(newman_pbmc$annotation$mixture))

colnames(newman_pbmc$annotation$mixture)[which(str_detect(colnames(newman_pbmc$annotation$mixture), "T cells regulatory"))] <- "T cells regulatory (Tregs)"
newman_pbmc$annotation$mixture <- newman_pbmc$annotation$mixture[colnames(LM22)]
cbind(colnames(LM22) , colnames(newman_pbmc$annotation$mixture))

##map between LM22 and TIL10
lm22.to.til10 <- c("B cells","B cells","B cells","T cells CD8","T cells CD4","T cells CD4","T cells CD4",
                   "T cells CD4","Tregs","T cells CD8","NK cells","NK cells","Monocytes","Macrophages M1",
                   "Macrophages M1","Macrophages M2","Dendritic cells","Dendritic cells","Other","Other","Other","Neutrophils") 
length(lm22.to.til10)
til10.ct.raw  <- colnames(TIL10)
til10.ct.raw <- str_replace_all(til10.ct.raw,"\\."," ")

PBMC.composite <- CompositeMixture(newman_pbmc$annotation$mixture[1:20,],lm22.to.til10)
PBMC.composite <- PBMC.composite[, til10.ct.raw]
PBMC.QUA <- data.matrix(quanti_PBMC[,-c(1,12)])
PBMC.mixture <- GetMixture(PBMC.robust,"prop")
PBMC.ABB <- GetMixture(PBMC.abbas,"prop")
PBMC.ABI <- GetMixture(abis_PBMC,"prop")
PBMC.cibersort <- PBMC.cib.site
TIL10.PBMC <- data.frame(betahat = c(c(PBMC.cibersort), c(PBMC.ABB), c(PBMC.mixture),
                                     c(PBMC.ABI), c(PBMC.QUA)) ,
                         Method = factor( rep(c("CIBERSORT", "ABBAS","MIXTURE","ABIS","QUANTISEQ"), 
                                              each = length(as.numeric(PBMC.composite))), 
                                          levels = c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")),
                         CT = rep(colnames(PBMC.composite),each= nrow(PBMC.composite),times = 5),
                         betasim = rep(as.numeric(PBMC.composite), times = 5))

TIL10.PBMC$difs <- TIL10.PBMC$betahat - TIL10.PBMC$betasim
TIL10.PBMC <- setbetacol(TIL10.PBMC)

ROCAnalisis(TIL10.PBMC)
#           TP TN FP FN Se Sp   PPV   NPV    DO DOSeSp  ERG       F1
# CIBERSORT 90 25 75 10 90 25 54.55 71.43 92.77  75.66 0.85 67.92805
# ABBAS     91 41 59  9 91 41 60.67 82.00 73.71  59.68 0.68 72.80240
# MIXTURE   86 33 67 14 86 33 56.21 70.21 86.54  68.45 0.81 67.98481
# ABIS      95 46 54  5 95 46 63.76 90.20 65.96  54.23 0.59 76.30637
# QUANTISEQ 27 94  6 73 27 94 81.82 56.29 87.21  73.25 0.79 40.60173

BA.PBMC <- BA.plot(TIL10.PBMC,plot=T)
Cor.PBMC <- Cor.plot(TIL10.PBMC,plot=T,scale = scale_y_continuous(
  breaks = c(-2,-1,0,1,2),
  label = c("-2.0","-1.0", "0.0", "1.0","2.0")))

ggplot(TIL10.PBMC,aes(x=difs, group=beta,fill=beta)) + geom_density(alpha=0.4,size = 0.2) + labs(x="error", fill = expression(beta))+facet_grid(~ Method) 

save(TIL10.FL, TIL10.PBMC, file = "flow_cyotmetry_output_TIL10.RData")
