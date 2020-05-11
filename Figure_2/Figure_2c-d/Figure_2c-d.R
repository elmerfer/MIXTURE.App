#Figure_2c-d: Performance against real flow cytometry data.
rm(list=ls())
library(ggplot2)
library(plyr)
library(gridExtra)
#Load and format Follicular Lymphoma and PBMC output using LM22 signature
load("Figure_2/Figure_2c-d/flow_cyotmetry_output_LM22.RData")

df.fl.lm22$Method<-gsub("QUANTISEQ","quanTIseq",df.fl.lm22$Method)
df.fl.lm22$Method<-factor(df.fl.lm22$Method,c(levels="ABBAS","ABIS","quanTIseq","CIBERSORT","MIXTURE"))
df.fl.lm22$CT<-gsub("T cells CD8         ","T cells CD8",gsub("T cells CD4                .","T cells CD4",gsub("B cells ","B cells",gsub("Other         ","Other",df.fl.lm22$CT))))
df.fl.lm22$CT<-factor(df.fl.lm22$CT,levels=c("T cells CD8","T cells CD4","B cells","Other"))
df.fl.lm22$Signature = "LM22"
df.fl.lm22$BD= "FL"

df.pbmc.lm22$Method<-gsub("QUANTISEQ","quanTIseq",df.pbmc.lm22$Method)
df.pbmc.lm22$Method<-factor(df.pbmc.lm22$Method,c(levels="ABBAS","ABIS","quanTIseq","CIBERSORT","MIXTURE"))
df.pbmc.lm22$Signature = "LM22"
df.pbmc.lm22$BD= "PBMC"
df.pbmc.lm22$CT<-factor(df.pbmc.lm22$CT,levels=c("B cells naive","B cells memory","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells gamma delta","Monocytes","NK cells activated","Other"))

#Load and format Follicular Lymphoma and PBMC output using TIL10 signature
load("Figure_2/Figure_2c-d/flow_cyotmetry_output_TIL10.RData")

df.fl.til10$Method<-gsub("QUANTISEQ","quanTIseq",df.fl.til10$Method)
df.fl.til10$Method<-factor(df.fl.til10$Method,c(levels="ABBAS","ABIS","quanTIseq","CIBERSORT","MIXTURE"))
df.fl.til10$Signature = "TIL10"
df.fl.til10$BD = "FL"

df.pbmc.til10$Method<-gsub("QUANTISEQ","quanTIseq",df.pbmc.til10$Method)
df.pbmc.til10$Method<-factor(df.pbmc.til10$Method,c(levels="ABBAS","ABIS","quanTIseq","CIBERSORT","MIXTURE"))
df.pbmc.til10$Signature = "TIL10"
df.pbmc.til10$BD="PBMC"
df.pbmc.til10$CT<-gsub("R","Other",df.pbmc.til10$CT)
df.pbmc.til10$CT<-factor(df.pbmc.til10$CT,levels=c("B cells","T cells CD8","T cells CD4","Monocytes","NK cells", "Other"))

#Define functions to generate plots
BA.plot.cytometry <- function(df, title = NULL, xlab = FALSE, plot=FALSE, scale = NULL, type = NULL){
  df.summary <- ddply(df, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
                      min= min(betahat,na.rm=T), max = max(betahat,na.rm=T), cor = cor(betahat, betasim))
  qq <-c(0.1,NA, 0.8*max(df$betasim))
  qq[2] <- qq[1]+ (qq[3]-qq[1])/2
  qq <- round(qq,2)
  if(type==1){
  p <- ggplot(df, aes(betasim, difs)) + 
    geom_point(size =0.7,aes(colour = CT)) + 
    geom_hline(yintercept = 0, colour = "blue",linetype = "dotted") +
    geom_hline(data=df.summary,aes(yintercept=round(mean,3)), color = "red", linetype = "solid") +
    geom_hline(data=df.summary,aes(yintercept=round(mean+2*sd,3)), color = "red", linetype = "dashed") + 
    geom_hline(data=df.summary,aes(yintercept=round(mean-2*sd,3)), color = "red", linetype = "dashed") +
    geom_smooth() +  scale_x_continuous(breaks = qq)+ theme_bw() + 
    theme(strip.text = element_text(colour = "black"),axis.ticks.x.bottom = element_line(colour = "black"),
          panel.grid = element_blank(),axis.title.y.left = element_text(size = 18),strip.text.x = element_text(size=18),
          axis.text.x = element_text(size=14,colour = "black"),axis.text.y = element_text(size=16),
          plot.title = element_text(hjust = 0.5,size = 20), legend.position = "bottom")+
    ylim(-.45, 0.6)
  
  if(!is.null(scale)) p <- p + scale
  if(!is.null(title)) p <- p + ggtitle(title)
  if(xlab) {
    p <- p + labs( x = "",y = paste("Error","\n"))
  }else p <- p + labs( x = "",y =  paste("Error","\n"))
  }else{
      p <- ggplot(df, aes(betasim, difs)) + 
        geom_point(size =0.7,aes(colour = CT)) + 
        geom_hline(yintercept = 0, colour = "blue",linetype = "dotted") +
        geom_hline(data=df.summary,aes(yintercept=round(mean,3)), color = "red", linetype = "solid") +
        geom_hline(data=df.summary,aes(yintercept=round(mean+2*sd,3)), color = "red", linetype = "dashed") + 
        geom_hline(data=df.summary,aes(yintercept=round(mean-2*sd,3)), color = "red", linetype = "dashed") +
        geom_smooth() +  scale_x_continuous(breaks = qq)+ theme_bw() + 
        theme(strip.text = element_text(colour = "black"),axis.ticks.x.bottom = element_line(colour = "black"),
              panel.grid = element_blank(),panel.grid.major = element_blank(),axis.title.y.left = element_text(size = 0),
              strip.text.x = element_text(size=18),axis.text.x = element_text(size=14,colour = "black"),
              axis.text.y = element_text(size=16),plot.title = element_text(hjust = 0.5,size=20), legend.position = "bottom")+
        ylim(-.45, 0.6)
      
      if(!is.null(scale)) p <- p + scale
      if(!is.null(title)) p <- p + ggtitle(title)
      if(xlab) {
        p <- p + labs( x = "Simulated coefficients.",y = "Error")
      }else p <- p + labs( x = "",y = "Error")
    }
  p <-  p +    facet_wrap(~Method, nrow = 1)
  if(plot) print(p)
  return(p)
}
Cor.plot.cytometry <- function(df, title = NULL, xlab = FALSE, plot=FALSE, scale=NULL, type=NULL){
  df.summary <- ddply(df, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
                      min= min(betahat,na.rm=T), max = max(betahat,na.rm=T), cor = cor(betahat, betasim))
  qq <-c(0.1,NA, 0.8*max(df$betasim))
  qq[2] <- qq[1]+ (qq[3]-qq[1])/2
  qq <- round(qq,2)
  if(type==1){
  p <- ggplot(df, aes(betasim, betahat)) + 
    geom_point(size =0.7, aes(colour = CT)) + 
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    geom_smooth(method = "lm") + scale_x_continuous(
      breaks = qq)+ theme_bw() + theme(strip.background.x = element_blank(),strip.text = element_text(colour = "white"),
                                       panel.grid = element_blank(),strip.text.x=element_text(size = 0),
                                       axis.title.y.left = element_text(size = 20),axis.text.x = element_text(size=16),
                                       axis.text.y = element_text(size=16), axis.title.x = element_text(size = 20),
                                       legend.position = "bottom",legend.text = element_text(size = 16),
                                       legend.title = element_text(size = 16))
  }else{
    p <- ggplot(df, aes(betasim, betahat)) + 
      geom_point(size =0.7, aes(colour = CT)) + 
      geom_abline(intercept = 0, slope = 1, colour = "red") +
      geom_smooth(method = "lm") + scale_x_continuous(
        breaks = qq)+ theme_bw() + theme(strip.background = element_blank(),strip.text = element_text(colour = "white"),
                                         strip.text.x=element_text(size = 0),axis.title.y = element_text(size = 0),
                                         axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), 
                                         legend.position = "bottom", axis.title.x = element_text(size = 20),
                                         panel.grid = element_blank(),legend.text = element_text(size = 16),
                                         legend.title = element_text(size = 16))
  }
  if(!is.null(scale)) p <- p +  scale
  if(!is.null(title)) p <- p + ggtitle(title)
  if(xlab) {
    p <- p + labs( x = "Flow Cytometry Coefficients", y = paste("Estimated","\ncoefficients"))
  }else p <- p + labs( x = "",y = paste("Estimated","\ncoefficients"))
  p <-  p +    facet_wrap(~Method, nrow = 1)
  if(plot) print(p)
  return(p)
}

#Bland-Altman plots for Follicular Lymphoma data
#LM22
BA.FL.LM22 <- BA.plot.cytometry(df.fl.lm22,plot=T,title = "LM22",type=1)
BA.FL.LM22[["labels"]][["colour"]]<-"Cell Type"
#TIL10
BA.FL.til10 <- BA.plot.cytometry(df.fl.til10,plot=T,title = "TIL10",type=2)
BA.FL.til10[["labels"]][["colour"]]<-"Cell Type"

#Correlation plots for Follicular Lymphoma data
#LM22
Cor.FL.LM22 <- Cor.plot.cytometry(df.fl.lm22,plot=T,type=1,scale = scale_y_continuous(breaks = c(-2,-1,0,1,2),label = c("-2.00","-1.00", "0.00", "1.00","2.00")))
Cor.FL.LM22[["labels"]][["colour"]]<-"Cell Type"
#TIL10
Cor.FL.til10 <- Cor.plot.cytometry(df.fl.til10,plot=T,type=2,scale = scale_y_continuous(breaks = c(-1,0,1),label = c("-1.00", "0.00", "1.00")))
Cor.FL.til10[["labels"]][["colour"]]<-"Cell Type"

#All plots for Follicular Lymphoma data
dev.off()
ggsave(filename = "Figure_2/Figure_2c-d/Figure_2c-d_FL.pdf",plot=grid.arrange(BA.FL.LM22,BA.FL.til10,Cor.FL.LM22,Cor.FL.til10,nrow=2,ncol=2),device = "pdf",width = unit(36,"cm"),height = unit(8,"cm"))

#Bland-Altman plots for PBMCs data
#LM22
BA.PBMC.LM22 <- BA.plot.cytometry(df.pbmc.lm22,plot=T,title="LM22",type=1)
BA.PBMC.LM22[["labels"]][["colour"]]<-"Cell Type"
#TIL10
BA.PBMC.til10 <- BA.plot.cytometry(df.pbmc.til10,plot=T,title="TIL10",type=2)
BA.PBMC.til10[["labels"]][["colour"]]<-"Cell Type"

#Correlation plots for PBMCs data
#LM22
Cor.PB.LM22 <- Cor.plot.cytometry(df.pbmc.lm22,plot=T,xlab = TRUE,type=1,scale = scale_y_continuous(breaks = c(-1,0,1),label = c("-1.00", "0.00", "1.00")))
Cor.PB.LM22[["labels"]][["colour"]]<-"Cell Type"
#TIL10
Cor.PB.til10 <- Cor.plot.cytometry(df.pbmc.til10,plot=T,xlab = TRUE,type=2,scale = scale_y_continuous(breaks = c(-1,0,1),label = c("-1.00", "0.00", "1.00")))
Cor.PB.til10[["labels"]][["colour"]]<-"Cell Type"

#Save all plots for Follicular Lymphoma data in PDF 36x8
dev.off()
ggsave(filename = "Figure_2/Figure_2c-d/Figure_2c-d_PBMC.pdf",plot=grid.arrange(BA.PBMC.LM22,BA.PBMC.til10,Cor.PB.LM22,Cor.PB.til10,nrow=2,ncol=2),device = "pdf",width = unit(36,"cm"),height = unit(8,"cm"))


