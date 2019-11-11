# Install immunedeconv
 # install.packages("remotes")
 # remotes::install_github("grst/immunedeconv")
rm(list=ls())
library(data.table)
library(ComplexHeatmap)

library(ggplot2)
library(circlize)
library(dtangle)
library(immunedeconv)
library(gridExtra)
source('Utils/MIXTURE.DEBUG_V0.1.R')

##get path to TIL10 signature file
load("Data/LM22.RData")
load("Data/TIL10.RData")

##Simulated data from LM22
betas.list <- readRDS("Data/betas.list.rds")
M.pure <- do.call(cbind, lapply(betas.list, function(x) x$A))
betas.noise.list <- readRDS("Data/betas.noise.list.rds")
M.pure.noise <- do.call(cbind, lapply(betas.noise.list, function(x) x$A)) 

#Simulated data from TIL10
TIL10.betas.list <- readRDS("Data/TIL10.betas.list.rds")
TIL10.M.pure <- do.call(cbind, lapply(TIL10.betas.list, function(x) x$A)) 
TIL10.betas.noise.list <- readRDS("Data/TIL10.betas.noise.list.rds")
TIL10.M.pure.noise <- do.call(cbind, lapply(TIL10.betas.noise.list, function(x) x$A)) 
 
 # Save data to run CIBERSORT web
 # aux <- data.frame("Gene Symbol" = rownames(TIL10.M.pure),TIL10.M.pure)
 # write.table(aux, file="/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/TIL10.M.pure.txt", quote = FALSE, row.names = FALSE, sep="\t")
 # aux <- data.frame("Gene Symbol" = rownames(TIL10.M.pure.noise),TIL10.M.pure.noise)
 # write.table(aux, file="/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/TIL10.M.pure.noise.txt", quote = FALSE, row.names = FALSE, sep="\t")
 
nc.cores <- 7L

 TIL10.pure.quanti <-deconvolute_quantiseq.default(mix.mat = TIL10.M.pure, 
                                              arrays = FALSE, 
                                              signame = "TIL10", 
                                              tumor = FALSE, 
                                              mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
 pure.quanti <-deconvolute_quantiseq.default(mix.mat = M.pure, 
                                                   arrays = FALSE, 
                                                   signame = "/home/elmer/Dropbox/IDEAS/cibersort/MIXTURE/Data/LM22", 
                                                   tumor = FALSE, 
                                                   mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
 
 TIL10.MIX <- MIXTURE(expressionMatrix = TIL10.M.pure, signatureMatrix =  TIL10, functionMixture =  nu.svm.robust.RFE, useCores = nc.cores)
 MIX <- MIXTURE(expressionMatrix = M.pure, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = nc.cores)
 
 
 TIL10.ABBAS.pure <- MIXTURE(expressionMatrix = TIL10.M.pure, signatureMatrix =  TIL10, functionMixture =  ls.rfe.abbas, useCores = nc.cores)
 ABBAS.pure <- MIXTURE(expressionMatrix = M.pure, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = nc.cores)
 
 # out.dt <- dtangle(Y = log2(t(LM22)) , reference = log2(t(LM22)))
 ##RLM
 TIL10.ABIS.pure <- MIXTURE(expressionMatrix = TIL10.M.pure, signatureMatrix =  TIL10, functionMixture =  rlm.abis, useCores = 1L)
 ABIS.pure <- MIXTURE(expressionMatrix = M.pure, signatureMatrix =  LM22, functionMixture =  rlm.abis, useCores = 1L)
 
 TIL10.noise.cib <- ReadCibersortWebResults(file="Data/CIBERSORT.Output_TIL10.Mix.Noise.csv", type = "csv", nct = 10)
 TIL10.pure.cib <- ReadCibersortWebResults(file="Data/CIBERSORT.Output_TIL10.Mix.pure.csv", type = "csv", nct = 10)
 
 M.r.cib <- ReadCibersortWebResults("Data/MixPure.CIBERSORT.Output_Job5.csv", type = "csv")
 M.r.cib.noise <-ReadCibersortWebResults("Data/NoiseMix.CIBERSORT.Output_Job4.csv", type = "csv")
 
 TIL10.pure.noise.quanti <-deconvolute_quantiseq.default(mix.mat = TIL10.M.pure.noise, 
                                                   arrays = FALSE, 
                                                   signame = "TIL10", 
                                                   tumor = FALSE, 
                                                   mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
 pure.noise.quanti <-deconvolute_quantiseq.default(mix.mat = M.pure.noise, 
                                             arrays = FALSE, 
                                             signame = "/home/elmer/Dropbox/IDEAS/cibersort/MIXTURE/Data/LM22", 
                                             tumor = FALSE, 
                                             mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
 
 TIL10.MIX.noise <- MIXTURE(expressionMatrix = TIL10.M.pure.noise, signatureMatrix =  TIL10, functionMixture =  nu.svm.robust.RFE, useCores = 10L)
 MIX.noise <- MIXTURE(expressionMatrix = M.pure.noise, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = 10L)
 
 TIL10.ABBAS.noise <- MIXTURE(expressionMatrix = TIL10.M.pure.noise, signatureMatrix =  TIL10, functionMixture =  ls.rfe.abbas, useCores = 3L)
 ABBAS.noise <- MIXTURE(expressionMatrix = M.pure.noise, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = 3L)
 
 # out.dt <- dtangle(Y = log2(t(LM22)) , reference = log2(t(LM22)))
 ##RLM
 TIL10.ABIS.noise <- MIXTURE(expressionMatrix = TIL10.M.pure.noise, signatureMatrix =  TIL10, functionMixture =  rlm.abis, useCores = 1L)
 ABIS.noise <- MIXTURE(expressionMatrix = M.pure.noise, signatureMatrix =  LM22, functionMixture =  rlm.abis, useCores = 1L)
 
 
 TIL10.Beta <- do.call(rbind, lapply(TIL10.betas.list, function(x) x$beta))
 TIL10.Beta.n <- do.call(rbind, lapply(TIL10.betas.noise.list, function(x) x$beta))
 Beta <- do.call(rbind, lapply(betas.list, function(x) x$beta))
 Beta.n <- do.call(rbind, lapply(betas.noise.list, function(x) x$beta))
 
 
 TIL10.df.BlandAltman <- data.frame(betahat = c( as.numeric(data.matrix(TIL10.pure.quanti[,-c(1,12)])), 
                                                            as.numeric(GetMixture(TIL10.MIX)),
                                                 as.numeric(GetMixture(TIL10.ABBAS.pure)),
                                                 as.numeric(GetMixture(TIL10.ABIS.pure)),
                                                 as.numeric(TIL10.pure.cib)),
                              betasim = rep(as.numeric(TIL10.Beta),5),
                              Method = factor(rep(c("QUANTISEQ","MIXTURE","ABBAS","ABIS","CIBERSORT"), 
                                                  each=length(as.numeric(TIL10.Beta))),
                                                  levels =  c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")))
 
 TIL10.df.BlandAltman$difs <- TIL10.df.BlandAltman$betahat - TIL10.df.BlandAltman$betasim
 
 df.BlandAltman <- data.frame(betahat = c( as.numeric(data.matrix(pure.quanti[,-c(1,24)])), 
                                                       as.numeric(GetMixture(MIX)),
                                                       as.numeric(GetMixture(ABBAS.pure)),
                                                       as.numeric(GetMixture(ABIS.pure)),
                                                       as.numeric(M.r.cib)),
                                          betasim = rep(as.numeric(Beta),5),
                                          Method = factor(rep(c("QUANTISEQ","MIXTURE","ABBAS","ABIS","CIBERSORT"), 
                                                              each=length(as.numeric(Beta))),
                                                          levels =  c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")))
 df.BlandAltman$difs <- df.BlandAltman$betahat - df.BlandAltman$betasim 
 
 TIL10.df.BlandAltman.n <- data.frame(betahat = c( as.numeric(data.matrix(TIL10.pure.noise.quanti[,-c(1,12)])), 
                                                 as.numeric(GetMixture(TIL10.MIX.noise)),
                                                 as.numeric(GetMixture(TIL10.ABBAS.noise)),
                                                 as.numeric(GetMixture(TIL10.ABIS.noise)),
                                                 as.numeric(TIL10.noise.cib)),
                                    betasim = rep(as.numeric(TIL10.Beta.n),5),
                                    Method = factor(rep(c("QUANTISEQ","MIXTURE","ABBAS","ABIS","CIBERSORT"), 
                                                        each=length(as.numeric(TIL10.Beta.n))),
                                    levels =  c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")))
 
 TIL10.df.BlandAltman.n$difs <- TIL10.df.BlandAltman.n$betahat - TIL10.df.BlandAltman.n$betasim
 
 
 df.BlandAltman.n <- data.frame(betahat = c( as.numeric(data.matrix(pure.noise.quanti[,-c(1,24)])), 
                                                   as.numeric(GetMixture(MIX.noise)),
                                                   as.numeric(GetMixture(ABBAS.noise)),
                                                   as.numeric(GetMixture(ABIS.noise)),
                                                   as.numeric(M.r.cib.noise)),
                                      betasim = rep(as.numeric(Beta.n),5),
                                      Method = factor(rep(c("QUANTISEQ","MIXTURE","ABBAS","ABIS","CIBERSORT"), 
                                                          each=length(as.numeric(Beta.n))),
                                                      levels =  c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")))
 df.BlandAltman.n$difs <- df.BlandAltman.n$betahat - df.BlandAltman.n$betasim 

  
 





dout <- ddply(TIL10.df.BlandAltman.n, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
      min= min(difs,na.rm=T), max = max(difs,na.rm=T), range = max(difs, na.rm=T)-min(difs, na.rm=T),cor = cor(betasim,betahat))
rownames(dout) <- dout[,1]
round(dout[,-1],3)
#           mean    sd    min   max range   cor
# ABBAS        0 0.083 -0.615 0.514 1.129 0.831
# ABIS         0 0.128 -0.377 0.718 1.095 0.667
# QUANTISEQ    0 0.266 -0.815 1.000 1.815 0.446
# CIBERSORT    0 0.078 -0.793 0.768 1.561 0.854
# MIXTURE      0 0.023 -0.243 0.341 0.585 0.988
dout <- ddply(subset(TIL10.df.BlandAltman.n, betasim >0), .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
              min= min(difs,na.rm=T), max = max(difs,na.rm=T), range = max(difs, na.rm=T)-min(difs, na.rm=T),cor = cor(betasim,betahat))
rownames(dout) <- dout[,1]
round(dout[,-1],3)
#           mean    sd    min   max range   cor
# ABBAS     -0.043 0.094 -0.615 0.487 1.103 0.751
# ABIS      -0.042 0.114 -0.377 0.572 0.949 0.672
# QUANTISEQ -0.019 0.399 -0.815 0.934 1.748 0.294
# CIBERSORT -0.029 0.088 -0.793 0.570 1.363 0.784
# MIXTURE   -0.008 0.026 -0.243 0.199 0.442 0.981
# 

dout <- ddply(subset(TIL10.df.BlandAltman.n, betasim <=0), .(Method), summarise, min= min(difs,na.rm=T), Q1= quantile(difs, 0.25),
              mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T),Q3= quantile(difs, 0.75),
              max = max(difs,na.rm=T), range = max(difs, na.rm=T)-min(difs, na.rm=T),cor = cor(betasim,betahat))
rownames(dout) <- dout[,1]
round(dout[,-1],3)
# mean    sd   min   max range cor
# ABBAS     0.029 0.059  0.00 0.514 0.514  NA
# ABIS      0.028 0.129 -0.21 0.718 0.928  NA
# QUANTISEQ 0.012 0.108  0.00 1.000 1.000  NA
# CIBERSORT 0.019 0.065  0.00 0.768 0.768  NA
# MIXTURE   0.005 0.020  0.00 0.341 0.341  NA

dout <- ddply(df.BlandAltman.n, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
              min= min(difs,na.rm=T), max = max(difs,na.rm=T), range = max(difs, na.rm=T)-min(difs, na.rm=T),cor = cor(betasim,betahat))
rownames(dout) <- dout[,1]
round(dout[,-1],3)
#           mean    sd    min   max range   cor
# ABBAS        0 0.046 -0.443 0.485 0.928 0.894
# ABIS         0 0.049 -0.202 0.480 0.682 0.882
# QUANTISEQ    0 0.063 -0.526 0.588 1.114 0.815
# CIBERSORT    0 0.011 -0.158 0.167 0.325 0.996
# MIXTURE      0 0.005 -0.113 0.107 0.220 0.999
dout <- ddply(subset(df.BlandAltman.n, betasim >0), .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
              min= min(difs,na.rm=T), max = max(difs,na.rm=T), range = max(difs, na.rm=T)-min(difs, na.rm=T),cor = cor(betasim,betahat))
rownames(dout) <- dout[,1]
round(dout[,-1],3)
#         mean    sd    min   max range   cor
# ABBAS     -0.039 0.060 -0.443 0.387 0.830 0.874
# ABIS      -0.022 0.048 -0.202 0.311 0.514 0.917
# QUANTISEQ -0.036 0.104 -0.526 0.495 1.021 0.755
# CIBERSORT -0.011 0.015 -0.158 0.167 0.325 0.994
# MIXTURE   -0.003 0.008 -0.113 0.052 0.164 0.998

BA.plot <- function(df, title = NULL, xlab = FALSE, plot=FALSE, scale = NULL){
  df.summary <- ddply(df, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
                      min= min(betahat,na.rm=T), max = max(betahat,na.rm=T), cor = cor(betahat, betasim))
  qq <-c(0.1,NA, 0.8*max(df$betasim))
  qq[2] <- qq[1]+ (qq[3]-qq[1])/2
  qq <- round(qq,2)
  p <- ggplot(df, aes(betasim, difs)) + 
    geom_point(size =0.5) + 
    geom_hline(yintercept = 0, colour = "blue",linetype = "dotted") +
    geom_hline(data=df.summary,aes(yintercept=round(mean,3)), color = "red", linetype = "solid") +
    geom_hline(data=df.summary,aes(yintercept=round(mean+2*sd,3)), color = "red", linetype = "dashed") + 
    geom_hline(data=df.summary,aes(yintercept=round(mean-2*sd,3)), color = "red", linetype = "dashed") +
    geom_smooth() +  scale_x_continuous(
      breaks = qq)+ theme_bw() + theme(axis.text.x = element_text(size=9),axis.text.y = element_text(size=9))+
    theme_bw()+  ylim(-.45, 0.6)
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
    geom_point(size =0.5) + 
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    geom_smooth(method = "lm") + scale_x_continuous(
      breaks = qq)+ theme_bw() + theme(axis.text.x = element_text(size=9),axis.text.y = element_text(size=9))
  
  if(!is.null(scale)) p <- p +  scale
  
  if(!is.null(title)) p <- p + ggtitle(title)
  if(xlab) {
    p <- p + labs( x = "Simulated coefficients.",y = "Est. coeffs")
  }else p <- p + labs( x = "",y = "Est. coeffs")
  p <-  p +    facet_wrap(~Method, nrow = 1)
  
  if(plot) print(p)
  return(p)
}

# setEPS()
# postscript("For Figure1.eps")
# grid.arrange(LM22.plot, TIL10.plot, nrow = 2)
# dev.off()

LM22.plot <- BA.plot(df.BlandAltman.n,"a")
LM22.cor <- Cor.plot(df.BlandAltman.n,"b",plot=F,scale = scale_y_continuous(
  breaks = c(0,0.5,1),
  label = c("0.00", "0.05","1.00")))

TIL10.plot <- BA.plot(TIL10.df.BlandAltman.n,"c")
TIL10.cor <- Cor.plot(TIL10.df.BlandAltman.n,title = "d",xlab = TRUE,scale = scale_y_continuous(
  breaks = c(0,0.5,1),
  label = c("0.00", "0.05","1.00")))

setEPS()
postscript("Correlation for Figure1.eps")
grid.arrange(LM22.cor, TIL10.cor, nrow = 2)
dev.off()
setEPS()
pdf("/home/elmer/Dropbox/IDEAS/cibersort/GenomeR/Figuras/FigurasFinales/FIG2.pdf", paper = "a4",
    width = 10, height = 15)
grid.arrange(LM22.plot, LM22.cor, TIL10.plot, TIL10.cor, nrow = 4)
dev.off()

pdf("FIG2.pdf",paper="a4", width = 10, height = 10)
grid.arrange(LM22.plot, LM22.cor, TIL10.plot, TIL10.cor, nrow = 4)
dev.off()



##Number of estimated cell types
beta.length.noise <- apply(Beta.noise ,1, function(x) sum(x>0))
TIL10.beta.length.noise <- apply(TIL10.Beta.n ,1, function(x) sum(x>0))

##number of estimated betas
cib.bet.length.noise <- rowSums(M.r.cib.noise>0)
TIL10.cib.bet.length.noise <- rowSums(TIL10.pure.cib>0)

robust.bet.length.noise2 <- rowSums(GetMixture(MIX.noise,"prop")>0)
TIL10.robust.bet.length.noise2 <- rowSums(GetMixture(TIL10.MIX.noise,"prop")>0)

abbas.bet.length.noise <- rowSums(GetMixture(ABBAS.noise,"prop")>0)
TIL10.abbas.bet.length.noise <- rowSums(GetMixture(TIL10.ABBAS.noise,"prop")>0)


# dt.length.noise <- rowSums(dt_out.noise$estimates > 0) 

abis.bet.length.noise <- rowSums(GetMixture(ABIS.noise,"prop")>0)
TIL10.abis.bet.length.noise <- rowSums(GetMixture(TIL10.ABIS.noise,"prop")>0)

quanti.bet.length.noise <- rowSums(pure.noise.quanti[,-c(1,24)] >0)
TIL10.quanti.bet.length.noise <- rowSums(TIL10.pure.noise.quanti[,-c(1,12)] >0)

df.betas.noise <- data.frame(est = c(cib.bet.length.noise,
                                     robust.bet.length.noise2,
                                     abbas.bet.length.noise,
                                     # dt.length.noise,
                                     abis.bet.length.noise,
                                     quanti.bet.length.noise),
                             beta = as.factor(c(beta.length.noise,beta.length.noise,beta.length.noise,beta.length.noise,beta.length.noise)),
                             method = factor(c(rep("CIBERSORT", length(cib.bet.length.noise)), 
                                               rep("MIXTURE",length(robust.bet.length.noise2)),
                                               rep("ABBAS",length(abbas.bet.length.noise)),
                                               # rep("DTANGLE",length(beta.length.noise)),
                                               rep("ABIS",length(abis.bet.length.noise)),
                                               rep("QUANTISEQ",length(quanti.bet.length.noise))
                             ),levels = c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")))
df.betas.noise$Signature <- "LM22"

TIL10.df.betas.noise <- data.frame(est = c(TIL10.cib.bet.length.noise,
                                           TIL10.robust.bet.length.noise2,
                                           TIL10.abbas.bet.length.noise,
                                     # dt.length.noise,
                                     TIL10.abis.bet.length.noise,
                                     TIL10.quanti.bet.length.noise),
                             beta = as.factor(rep(TIL10.beta.length.noise,5)),
                             method = factor(c(rep("CIBERSORT", length(TIL10.cib.bet.length.noise)), 
                                               rep("MIXTURE",length(TIL10.robust.bet.length.noise2)),
                                               rep("ABBAS",length(TIL10.abbas.bet.length.noise)),
                                               # rep("DTANGLE",length(beta.length.noise)),
                                               rep("ABIS",length(TIL10.abis.bet.length.noise)),
                                               rep("QUANTISEQ",length(TIL10.quanti.bet.length.noise))
                             ),levels = c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")))
TIL10.df.betas.noise$Signature <- "TIL10"





bx1 <- ggplot(df.betas.noise, aes(x=beta, y=est)) +
  geom_boxplot(position = position_dodge(width = 0.7), aes(colour = method) ) +
  labs( x = "True number of coefficients",
        y = "Estimated number of coefficients", 
        fill = "Algorithm") + theme_bw() +
  scale_y_continuous(name ="Estimated number of coefficients", 
                    breaks = c(2:15),labels=as.character(c(2:15)))#+ facet_wrap(~Signature) 
scale_y_continuous(name, breaks, labels, limits, trans)
bx2 <- ggplot(TIL10.df.betas.noise, aes(x=beta, y=est)) +
  geom_boxplot(position = position_dodge(width = 0.7), aes(colour = method) ) +
  labs( x = "True number of coefficients",
        y = "Estimated number of coefficients", 
        fill = "Algorithm") + theme_bw() +
  scale_y_discrete(name ="Estimated number of coefficients", 
                   limits=as.character(c(2:10)))#+ facet_wrap(~Signature) 
setEPS()
postscript("BoxplotsEstimationNumber.eps")
grid.arrange(bx1, bx2, nrow = 2)
dev.off()


# save(df.betas.noise, file = "Data/DF.Betas.Noise.LM22.RData")
# save(TIL10.df.betas.noise, file = "Data/DF.Betas.Noise.TIL10.RData")



