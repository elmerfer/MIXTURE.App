#Figure_2c-d and Supplementary Figure_1: Performance against real flow cytometry data.
library(plyr)
library(gridExtra)
library(ggplot2)
library(ggpubr)

rm(list=ls())

#Load and generate usefull functions
source("Utils/dm_utils.R")
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
BA.plot <- function(df, title = NULL, xlab = FALSE, plot = FALSE, scale = NULL,legend = FALSE){
  df.summary <- ddply(df, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
                      min= min(betahat,na.rm=T), max = max(betahat,na.rm=T), cor = cor(betahat, betasim))
  qq <-c(0.1,NA, 0.8*max(df$betasim))
  qq[2] <- qq[1]+ (qq[3]-qq[1])/2
  qq <- round(qq,2)
  p <- ggplot(df, aes(betasim, difs)) + 
    geom_point(size =0.7,aes(colour = CT),show.legend = legend) + labs(colour="Cell Type") +
    geom_hline(yintercept = 0, colour = "blue",linetype = "dotted") +
    geom_hline(data=df.summary,aes(yintercept=round(mean,3)), color = "red", linetype = "solid") +
    geom_hline(data=df.summary,aes(yintercept=round(mean+2*sd,3)), color = "red", linetype = "dashed") + 
    geom_hline(data=df.summary,aes(yintercept=round(mean-2*sd,3)), color = "red", linetype = "dashed") +
    geom_smooth(aes(betasim, difs), data=subset(df,betasim>0)) +  scale_x_continuous(breaks = qq) +
    theme_bw() + theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),
                       legend.text = element_text(size = 16),axis.text.x = element_text(size=16),
                       plot.title = element_text(hjust = 0.5,size = 30),axis.text.y = element_text(size=16),
                       strip.text = element_text(size = 26),panel.grid = element_blank(),
                       legend.title = element_text(size=16,face="bold"),legend.position = "bottom") +
    ylim(-.45, 0.6)
  
  if(!is.null(scale)) {p <- p + scale}
  if(!is.null(title)) {p <- p + ggtitle(title)}
  if(xlab) {
    p <- p + labs( x = "Flow Cytometry Coefficients",y = "Error")
    }else{
      p <- p + labs( x = "",y = "Error")
    }
  p <-  p + facet_wrap(~Method,nrow = 1)
  if(plot) {print(p)}
  return(p)
}
Cor.plot <- function(df, title = NULL, xlab = FALSE, plot = FALSE, scale = NULL,legend = FALSE){
  df.summary <- ddply(df, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
                      min= min(betahat,na.rm=T), max = max(betahat,na.rm=T), cor = cor(betahat, betasim))
  qq <-c(0.1,NA, 0.8*max(df$betasim))
  qq[2] <- qq[1]+ (qq[3]-qq[1])/2
  qq <- round(qq,2)
  
  p <- ggplot(df, aes(betasim, betahat)) + 
    geom_point(size =0.7, aes(colour = CT),show.legend = legend) + labs(colour="Cell Type") +
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    geom_smooth(method = "lm") + scale_x_continuous(breaks = qq) + 
    theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust = 0.5,size = 30),
                       axis.title.y = element_text(size = 20),axis.text.x = element_text(size=16),
                       axis.text.y = element_text(size=16), axis.title.x = element_text(size = 20),
                       legend.position = "bottom",legend.text = element_text(size = 16),
                       strip.text = element_text(size = 26),legend.title = element_text(size = 16,face = "bold"))
  
  if(!is.null(scale)) {p <- p + scale}
  if(!is.null(title)) {p <- p + ggtitle(title)}
  if(xlab) {
    p <- p + labs( x = "Flow Cytometry Coefficients",y = "Estimated Coefficients")
  }else{
    p <- p + labs( x = "",y = "Estimated Coefficients")
  }
  p <-  p + facet_wrap(~Method,nrow = 1)
  if(plot) {print(p)}
  return(p)
}
Summary.performace <- function(df){
  cat("\n Betas == 0 \n")
  df.sum.ceros <- ddply(subset(df, betasim == 0), .(Method), summarise, min= min(difs, na.rm=T), q2 = quantile(difs, 0.25),
                        mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), med=median(difs, na.rm=T), max = max(difs, na.rm=T), min= min(difs, na.rm=T), q3 = quantile(difs, 0.75),
                        cor = cor(betasim,betahat))
  rownames(df.sum.ceros) <- df.sum.ceros[,1]
  print(round(df.sum.ceros[,-1],3)  )
  cat("\n Betas > 0 \n")
  df.sum.ceros <- ddply(subset(df, betasim > 0), .(Method), summarise, min= min(difs, na.rm=T), q2 = quantile(difs, 0.25),
                        mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), med=median(difs, na.rm=T), max = max(difs, na.rm=T), min= min(difs, na.rm=T), q3 = quantile(difs, 0.75),
                        cor = cor(betasim,betahat))
  rownames(df.sum.ceros) <- df.sum.ceros[,1]
  print(round(df.sum.ceros[,-1],3)  )
  cat("\n Total \n")
  df.sum.ceros <- ddply(df, .(Method), summarise, min= min(difs, na.rm=T), q2 = quantile(difs, 0.25),
                        mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), med=median(difs, na.rm=T), max = max(difs, na.rm=T), min= min(difs, na.rm=T), q3 = quantile(difs, 0.75),
                        cor = cor(betasim,betahat))
  rownames(df.sum.ceros) <- df.sum.ceros[,1]
  print(round(df.sum.ceros[,-1],3)  )
}

#Load Follicular Lymphoma and PBMC's data (Immune infiltrate estimation and Flow Cytometry data)
load(file = "Figure_2/Figure_2c-d_sup1/flow_cyotmetry_output_LM22.RData")
load(file = "Figure_2/Figure_2c-d_sup1/flow_cyotmetry_output_TIL10.RData")

###Working with LM22 signature

##Follicular Lymphoma dataset

#Bland-Altman plot
LM22.BA.FL <- BA.plot(df = LM22.FL,title="LM22",xlab=T,plot=T)

#Correlation plot
LM22.Cor.FL <- Cor.plot(LM22.FL,title="LM22",plot=T,xlab=T,scale = scale_y_continuous(
  breaks = c(-1,0,1),
  label = c("-1,00","0","1.00")))

#Error rates
error_rates_FL_LM22<- ggplot(LM22.FL,aes(x=difs, group=beta,fill=beta)) + geom_density(alpha=0.4,size = 0.2) + 
  labs(x="Error",fill = expression(beta),y="Density") + facet_grid(~ Method) + ggtitle("Error Rates - FL with LM22") +
  theme(strip.text = element_text(colour = "black"),strip.text.x=element_text(size = 14),
        axis.title.y.left = element_text(size = 14),axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12), axis.title.x = element_text(size = 14),
        legend.position = "right",legend.text = element_text(size = 14),
        legend.title = element_text(size = 14,face = "bold"))

##PBMC dataset

#Bland-Altman plot
LM22.BA.PBMC <- BA.plot(df = LM22.PBMC,title="LM22",xlab=T,plot=T)

#Correlation plot
LM22.Cor.PBMC <- Cor.plot(LM22.PBMC,plot=T,title="LM22",xlab=T,scale = scale_y_continuous(
  breaks = c(-1,0,0.75),
  label = c("-1,00","0","0.75")))

#Error rates
error_rates_PBMC_LM22 <- ggplot(LM22.PBMC,aes(x=difs, group=beta,fill=beta)) + geom_density(alpha=0.4,size = 0.2) + 
  labs(x="Error",fill = expression(beta),y="Density") + facet_grid(~ Method) + ggtitle("Error Rates - FL with LM22") +
  theme(strip.text = element_text(colour = "black"),strip.text.x=element_text(size = 14),
        axis.title.y.left = element_text(size = 14),axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12), axis.title.x = element_text(size = 14),
        legend.position = "right",legend.text = element_text(size = 14),
        legend.title = element_text(size = 14,face = "bold"))

###Working with LM22 signature

##Follicular Lymphoma dataset

#Bland-Altman plot
TIL10.BA.FL <- BA.plot(TIL10.FL,title="TIL10",xlab=T,plot=T)

#Correlation plot
TIL10.Cor.FL <- Cor.plot(TIL10.FL,title="TIL10",xlab=T,plot=T,scale = scale_y_continuous(
  breaks = c(-2,-1,0,1,2),
  label = c("-2.00","-1,00","0","1.00","2.00")))

#Error rates
error_rates_FL_TIL10 <- ggplot(TIL10.FL,aes(x=difs, group=beta,fill=beta)) + geom_density(alpha=0.4,size = 0.2) + 
  labs(x="Error",fill = expression(beta),y="Density") + facet_grid(~ Method) + ggtitle("Error Rates - FL with LM22") +
  theme(strip.text = element_text(colour = "black"),strip.text.x=element_text(size = 14),
        axis.title.y.left = element_text(size = 14),axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12), axis.title.x = element_text(size = 14),
        legend.position = "right",legend.text = element_text(size = 14),
        legend.title = element_text(size = 14,face = "bold"))

##PBMC dataset

#Bland-Altman plot
TIL10.BA.PBMC <- BA.plot(TIL10.PBMC,title="TIL10",xlab=T,plot=T)

#Correlation plot
TIL10.Cor.PBMC <- Cor.plot(TIL10.PBMC,title="TIL10",xlab=T,plot=T,scale = scale_y_continuous(
  breaks = c(-1,0,1),
  label = c("-1,00","0","1.00")))

#Error rates
error_rates_PBMC_TIL10 <- ggplot(TIL10.PBMC,aes(x=difs, group=beta,fill=beta)) + geom_density(alpha=0.4,size = 0.2) + 
  labs(x="Error",fill = expression(beta),y="Density") + facet_grid(~ Method) + ggtitle("Error Rates - FL with LM22") +
  theme(strip.text = element_text(colour = "black"),strip.text.x=element_text(size = 14),
        axis.title.y.left = element_text(size = 14),axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12), axis.title.x = element_text(size = 14),
        legend.position = "right",legend.text = element_text(size = 14),
        legend.title = element_text(size = 14,face = "bold"))

#Save Figure_2c-d in pdf 20x35(cm)
dev.off()
ggsave(filename = "Figure_2/Figure_2c-d_sup1/Figure2_c-d.pdf",plot=grid.arrange(LM22.BA.FL,TIL10.BA.FL,LM22.Cor.FL,TIL10.Cor.FL,LM22.BA.PBMC,TIL10.BA.PBMC,LM22.Cor.PBMC,TIL10.Cor.PBMC,nrow=4,ncol=2),device = "pdf",height = unit(20,"cm"),width = unit(35,"cm"))

#Save Supplementary Figure 1
dev.off()
ggsave(filename = "Figure_2/Figure_2c-d_sup1/Supplementary_Figure1.pdf",plot=grid.arrange(error_rates_FL_LM22,error_rates_FL_TIL10,error_rates_PBMC_LM22,error_rates_PBMC_TIL10,ncol=1,nrow=4),device = "pdf",width = unit(10,"cm"),height = unit(20,"cm"))


## Summaries
Summary.performace(LM22.FL)
# Betas == 0 
#              min    q2  mean    sd    med   max    q3 cor
# ABBAS      0.000  0.00 0.042 0.076  0.000 0.297 0.035  NA
# ABIS      -0.061 -0.01 0.032 0.072 -0.004 0.241 0.060  NA
# QUANTISEQ  0.000  0.00 0.011 0.039  0.000 0.259 0.000  NA
# CIBERSORT  0.000  0.00 0.025 0.038  0.008 0.205 0.038  NA
# MIXTURE    0.000  0.00 0.018 0.042  0.000 0.250 0.010  NA
# 
# Betas > 0 
#             min     q2   mean    sd    med   max     q3   cor
# ABBAS     -0.427 -0.214 -0.125 0.121 -0.096 0.103 -0.036 0.941
# ABIS      -0.403 -0.180 -0.096 0.134 -0.079 0.131  0.012 0.920
# QUANTISEQ -0.333 -0.123 -0.034 0.147 -0.031 0.384  0.052 0.949
# CIBERSORT -0.442 -0.195 -0.074 0.144 -0.007 0.090  0.031 0.961
# MIXTURE   -0.358 -0.116 -0.053 0.124 -0.010 0.114  0.033 0.957
# 
# Total 
#             min     q2 mean    sd    med   max    q3   cor
# ABBAS     -0.427  0.000    0 0.115  0.000 0.297 0.013 0.861
# ABIS      -0.403 -0.022    0 0.106 -0.005 0.241 0.032 0.879
# QUANTISEQ -0.333  0.000    0 0.083  0.000 0.384 0.000 0.943
# CIBERSORT -0.442  0.000    0 0.090  0.007 0.205 0.034 0.954
# MIXTURE   -0.358  0.000    0 0.078  0.000 0.250 0.019 0.953
Summary.performace(LM22.PBMC)
# Betas == 0 
#             min     q2  mean    sd    med   max    q3 cor
# ABBAS      0.000  0.000 0.033 0.062  0.000 0.304 0.042  NA
# ABIS      -0.091 -0.017 0.013 0.066 -0.005 0.323 0.020  NA
# QUANTISEQ  0.000  0.000 0.022 0.078  0.000 0.520 0.000  NA
# CIBERSORT  0.000  0.000 0.014 0.029  0.000 0.267 0.018  NA
# MIXTURE    0.000  0.000 0.007 0.026  0.000 0.290 0.000  NA
# 
# Betas > 0 
#             min     q2   mean    sd    med   max     q3   cor
# ABBAS     -0.396 -0.105 -0.047 0.120 -0.030 0.281  0.000 0.274
# ABIS      -0.411 -0.103 -0.019 0.151 -0.027 0.404  0.044 0.358
# QUANTISEQ -0.396 -0.109 -0.032 0.168 -0.041 0.546 -0.009 0.175
# CIBERSORT -0.326 -0.060 -0.021 0.093 -0.009 0.266  0.022 0.573
# MIXTURE   -0.326 -0.052 -0.010 0.098 -0.008 0.271  0.032 0.624
# 
# Total 
# min     q2 mean    sd    med   max    q3   cor
# ABBAS     -0.396 -0.017    0 0.098  0.000 0.304 0.029 0.288
# ABIS      -0.411 -0.041    0 0.110 -0.008 0.404 0.026 0.442
# QUANTISEQ -0.396 -0.022    0 0.125  0.000 0.546 0.000 0.265
# CIBERSORT -0.326  0.000    0 0.066  0.000 0.267 0.019 0.672
# MIXTURE   -0.326  0.000    0 0.066  0.000 0.290 0.002 0.723
Summary.performace(TIL10.FL)
# Betas == 0 
#             min     q2   mean    sd    med   max     q3 cor
# ABBAS      0.000  0.000  0.020 0.054  0.000 0.287  0.000  NA
# ABIS      -1.826 -0.305 -0.078 0.627 -0.198 3.314 -0.032  NA
# QUANTISEQ  0.000  0.000  0.000 0.000  0.000 0.000  0.000  NA
# CIBERSORT  0.000  0.000  0.054 0.089  0.000 0.348  0.083  NA
# MIXTURE    0.000  0.000  0.029 0.071  0.000 0.312  0.000  NA
# 
# Betas > 0 
#             min     q2   mean    sd    med   max     q3   cor
# ABBAS     -0.320 -0.117 -0.033 0.134 -0.029 0.212  0.053 0.932
# ABIS      -0.328 -0.091  0.131 0.320 -0.038 1.072  0.333 0.695
# QUANTISEQ -0.402 -0.089  0.000 0.188 -0.029 0.460  0.116 0.949
# CIBERSORT -0.384 -0.133 -0.090 0.118 -0.037 0.076 -0.014 0.951
# MIXTURE   -0.290 -0.077 -0.048 0.111 -0.021 0.198  0.011 0.946
# 
# Total 
#             min     q2 mean    sd    med   max    q3   cor
# ABBAS     -0.320 -0.010    0 0.096  0.000 0.287 0.000 0.939
# ABIS      -1.826 -0.243    0 0.541 -0.125 3.314 0.299 0.498
# QUANTISEQ -0.402 -0.003    0 0.114  0.000 0.460 0.000 0.953
# CIBERSORT -0.384 -0.026    0 0.122  0.000 0.348 0.024 0.895
# MIXTURE   -0.290 -0.002    0 0.095  0.000 0.312 0.000 0.933

Summary.performace(TIL10.PBMC)
# Betas == 0 
#             min     q2  mean    sd   med   max    q3 cor
# ABBAS      0.00  0.000 0.052 0.062 0.030 0.230 0.092  NA
# ABIS      -0.12 -0.032 0.009 0.057 0.003 0.185 0.042  NA
# QUANTISEQ  0.00  0.000 0.022 0.115 0.000 1.000 0.000  NA
# CIBERSORT  0.00  0.001 0.055 0.064 0.037 0.282 0.084  NA
# MIXTURE    0.00  0.000 0.055 0.060 0.041 0.252 0.085  NA
# 
# Betas > 0 
#             min     q2   mean    sd    med   max     q3    cor
# ABBAS     -0.341 -0.109 -0.052 0.121 -0.050 0.420 -0.012  0.645
# ABIS      -0.366 -0.102 -0.009 0.182 -0.045 0.641  0.054  0.598
# QUANTISEQ -0.620 -0.211 -0.022 0.385 -0.121 0.955 -0.036 -0.138
# CIBERSORT -0.440 -0.132 -0.055 0.169 -0.068 0.568  0.000  0.419
# MIXTURE   -0.431 -0.131 -0.055 0.148 -0.057 0.353  0.018  0.540
# 
# Total 
#           min     q2    mean    sd    med   max    q3   cor
# ABBAS     -0.341 -0.050    0 0.109  0.000 0.420 0.054 0.669
# ABIS      -0.366 -0.058    0 0.135 -0.015 0.641 0.045 0.700
# QUANTISEQ -0.620 -0.118    0 0.284  0.000 1.000 0.000 0.118
# CIBERSORT -0.440 -0.067    0 0.139  0.000 0.568 0.067 0.499
# MIXTURE   -0.431 -0.057    0 0.125  0.000 0.353 0.066 0.585
