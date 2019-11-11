library(ggplot2)
load("FL_PBMC_LM22.RData")
load("FL-PBMC_TIL10.RData")


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
    geom_smooth() +  scale_x_continuous(breaks = qq)+ theme_bw() + 
    theme(axis.text.x = element_text(size=9),axis.text.y = element_text(size=9), legend.position = "none")+
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


df.fl.lm22$Signature = "LM22"
df.fl.lm22$BD= "FL"

df.fl.til10$Signature = "TIL10"
df.fl.til10$BD = "FL"

df.pbmc.lm22$Signature = "LM22"
df.pbmc.lm22$BD= "PBMC"
df.pbmc.til10$Signature = "TIL10"
df.pbmc.til10$BD="PBMC"
df.total <- rbind(df.fl.lm22, df.fl.til10, df.pbmc.lm22, df.pbmc.til10)

pl <- ggplot(df.total, aes(Method, difs, fill=Method)) + 
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) + geom_hline(yintercept = 0, linetype="dotted", col="red")+
  theme_bw() + labs( x = "Methods",y = "Error")+ theme(axis.text.x = element_text(size=9, angle = 45)) +
  facet_grid(rows=vars(Signature), cols= vars(BD))

BA.FL <- BA.plot(df.fl.lm22,plot=T)
BA.FL.til10 <- BA.plot(df.fl.til10,plot=T)
Cor.FL <- Cor.plot(df.fl.lm22,plot=T,scale = scale_y_continuous(
  breaks = c(-2,-1,0,1,2),
  label = c("-2.00","-1.00", "0.00", "1.00","2.00")))
Cor.FL.til10 <- Cor.plot(df.fl.til10,plot=T,scale = scale_y_continuous(
  breaks = c(-2,-1,0,1,2),
  label = c("-2.00","-1.00", "0.00", "1.00","2.00")))

BA.PB <- BA.plot(df.pbmc.lm22,plot=T)
BA.PB.til10 <- BA.plot(df.pbmc.til10,plot=T)
Cor.PB <- Cor.plot(df.pbmc.lm22,plot=T,scale = scale_y_continuous(
  breaks = c(-2,-1,0,1,2),
  label = c("-2.00","-1.00", "0.00", "1.00","2.00")))
Cor.PB.til10 <- Cor.plot(df.pbmc.til10,plot=T,scale = scale_y_continuous(
  breaks = c(-2,-1,0,1,2),
  label = c("-2.00","-1.00", "0.00", "1.00","2.00")))

pdf("/home/elmer/Dropbox/IDEAS/cibersort/GenomeR/Figuras/LM22_FL_PBMC.pdf",paper = "a4", height = 10)
  grid.arrange(BA.FL, Cor.FL,BA.PB, Cor.PB, nrow=4)
  dev.off()
  pdf("/home/elmer/Dropbox/IDEAS/cibersort/GenomeR/Figuras/TIL10_FL_PBMC.pdf",paper = "a4", height = 10)
  grid.arrange(BA.FL.til10, Cor.FL.til10,BA.PB.til10, Cor.PB.til10, nrow=4)
  dev.off()

table(df.pbmc.lm22$CT)
library(ggplot2)
setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/GenomeR/Figuras/BA_Errors.eps")
print(pl)
dev.off()

pdf("/home/elmer/Dropbox/IDEAS/cibersort/GenomeR/Figuras/BA_Errors.pdf", paper = "a4", width = 10, height = 10)
print(pl)
dev.off()
