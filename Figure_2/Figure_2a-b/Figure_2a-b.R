#Figure_2a-b: Bland-Altman and correlation plots for simulated data.
rm(list=ls())
library(ggplot2)
library(plyr)
library(gridExtra)

#Load matrix with the alogirthms output using LM22 signature on simulated data with noise
LM22.BlandAltman.n <- readRDS(file = "Figure_2/Figure_2a-b/DF_BlandAtman_Noise_LM22.RDS")
LM22.BlandAltman.n$Method <- gsub("QUANTISEQ","quanTIseq",LM22.BlandAltman.n$Method)
LM22.BlandAltman.n$Method <- factor(LM22.BlandAltman.n$Method,levels=c("ABBAS","ABIS","quanTIseq","CIBERSORT","MIXTURE"))

#Load matrix with the alogirthms output using TIl10 signature on simulated data with noise
#Load matrix with MIXTURE output using TIL10 signature on simulated data
TIl10.BlandAltman.n <- readRDS(file = "Figure_2/Figure_2a-b/DF_BlandAltman_Noise_TIL10.RDS")
TIl10.BlandAltman.n$Method <- gsub("QUANTISEQ","quanTIseq",TIl10.BlandAltman.n$Method)
TIl10.BlandAltman.n$Method <- factor(TIl10.BlandAltman.n$Method,levels=c("ABBAS","ABIS","quanTIseq","CIBERSORT","MIXTURE"))

#Define functions to plot different subsets of data
BA.plot <- function(df, title = NULL, xlab = FALSE, plot=FALSE, scale = NULL, graph = NULL){
  df.summary <- ddply(df, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
                      min= min(betahat,na.rm=T), max = max(betahat,na.rm=T), cor = cor(betahat, betasim))
  qq <-c(0.1,NA, 0.8*max(df$betasim))
  qq[2] <- qq[1]+ (qq[3]-qq[1])/2
  qq <- round(qq,2)
  if(graph == "LM22"){
  p <- ggplot(df, aes(betasim, difs)) + 
    geom_point(size =0.5) + 
    geom_hline(yintercept = 0, colour = "blue",linetype = "dotted") +
    geom_hline(data=df.summary,aes(yintercept=round(mean,3)), color = "red", linetype = "solid") +
    geom_hline(data=df.summary,aes(yintercept=round(mean+2*sd,3)), color = "red", linetype = "dashed") + 
    geom_hline(data=df.summary,aes(yintercept=round(mean-2*sd,3)), color = "red", linetype = "dashed") +
    geom_smooth() +  scale_x_continuous(breaks = qq)+ 
    theme_bw() + theme(axis.title.x = element_text(size=22,colour = "white"),
                       axis.title.y = element_text(size=22),
                       axis.text.x.bottom = element_text(size=22),
                       axis.text.y.left = element_text(size=22),
                       plot.title = element_text(hjust = 0.5,size = 30),
                       axis.title = element_text(size = 30),
                       strip.text = element_text(size = 26),
                       panel.grid = element_blank())+
    ylim(-.45, 0.6)
  }else{
    p <- ggplot(df, aes(betasim, difs)) + 
      geom_point(size =0.5) + 
      geom_hline(yintercept = 0, colour = "blue",linetype = "dotted") +
      geom_hline(data=df.summary,aes(yintercept=round(mean,3)), color = "red", linetype = "solid") +
      geom_hline(data=df.summary,aes(yintercept=round(mean+2*sd,3)), color = "red", linetype = "dashed") + 
      geom_hline(data=df.summary,aes(yintercept=round(mean-2*sd,3)), color = "red", linetype = "dashed") +
      geom_smooth() +  scale_x_continuous(breaks = qq)+ 
      theme_bw() + theme(axis.title.x = element_text(size=22,colour = "white"),
                         axis.title.y = element_text(size=0),
                         axis.text.x.bottom = element_text(size=22),
                         axis.text.y.left = element_text(size=22),
                         plot.title = element_text(hjust = 0.5,size = 30),
                         axis.title = element_text(size = 30),
                         strip.text = element_text(size = 26),
                         panel.grid = element_blank())+
      ylim(-.45, 0.6)}
  p + ggtitle(title)+
  if(!is.null(scale)) p <- p + scale
  if(!is.null(title)) p <- p + ggtitle(title)
  if(xlab) {
    p <- p + labs( x = "Simulated coefficients.",y = "Error")
  }else p <- p + labs( x = "",y = "Error")
  p <-  p +    facet_wrap(~Method, nrow = 1)
  if(plot) print(p)
  return(p)
}
Cor.plot <- function(df, title = NULL, xlab = FALSE, plot=FALSE, scale=NULL, graph = NULL){
  df.summary <- ddply(df, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
                      min= min(betahat,na.rm=T), max = max(betahat,na.rm=T), cor = cor(betahat, betasim))
  qq <-c(0.1,NA, 0.8*max(df$betasim))
  qq[2] <- qq[1]+ (qq[3]-qq[1])/2
  qq <- round(qq,2)
  if(graph == "LM22"){
  p <- ggplot(df, aes(betasim, betahat)) + 
    geom_point(size =0.5) + 
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    geom_smooth(method = "lm") + scale_x_continuous(breaks = qq)+ 
    theme_bw() + theme(panel.grid = element_blank(),axis.title.x = element_text(size=22),axis.title.y = element_text(size=22),
                       axis.text.x.bottom = element_text(size=22),axis.title = element_text(size = 30),
                       axis.text.y.left = element_text(size=22),strip.text = element_text(size=26),
                       plot.title = element_text(hjust = 0.5,size = 30,colour = "white"))
  }else{
  p<- ggplot(df, aes(betasim, betahat)) + 
  geom_point(size =0.5) + 
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  geom_smooth(method = "lm") + scale_x_continuous(breaks = qq)+ theme_bw() + 
  theme(panel.grid = element_blank(),axis.title.x = element_text(size=22),axis.title.y = element_text(size=0),
        axis.text.x.bottom = element_text(size=22),axis.title = element_text(size = 30),axis.text.y.left = element_text(size=22),
        strip.text = element_text(size=26),plot.title = element_text(hjust = 0.5,size = 30,colour = "white"))}
  if(!is.null(scale)) p <- p +  scale
  if(!is.null(title)) p <- p + ggtitle(title)
  if(xlab) {
    p <- p + labs( x = "Simulated coefficients.",y = "Estimated coefficients")
  }else p <- p + labs( x = "",y = "Est. coeffs")
  p <-  p +    facet_wrap(~Method, nrow = 1)
  if(plot) print(p)
  return(p)
}

#Generate Bland-Altamn and correlation plots for LM22 signature
LM22_BA <- BA.plot(LM22.BlandAltman.n,title = "LM22",graph = "LM22")
LM22_cor <- Cor.plot(LM22.BlandAltman.n,xlab = TRUE,title = "LM22",scale = scale_y_continuous(breaks = c(0,0.5,1),label = c("0.00", "0.50","1.00")), graph = "LM22")

#Generate Bland-Altamn and correlation plots for TIL10 signature
TIL10_BA <- BA.plot(TIl10.BlandAltman.n,title="TIL10", graph = "TIL10")
TIL10_cor <- Cor.plot(TIl10.BlandAltman.n,title = "TIL10",xlab = TRUE,scale = scale_y_continuous(breaks = c(0,0.5,1),label = c("0.00", "0.50","1.00")),graph = "TIL10")

#Save plot in PDF 36x8
ggsave(filename = "Figure_2/Figure_2a-b/Figure_2a-b_simulated.pdf",plot=grid.arrange(LM22_BA,TIL10_BA,LM22_cor,TIL10_cor, nrow = 2,ncol=2),device = "pdf",width = unit(36,"cm"),height = unit(8,"cm"))

