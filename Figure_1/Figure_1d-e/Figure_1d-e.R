#Figure_1d-e: Number of cell types estimated by each algorithm
rm(list=ls())
library(dplyr)
library(gridExtra)
library(ggplot2)

#Load and format matrix with the output for each algorithm using LM22 signature of real and estimated number of cell types
LM22_betas_noise <- readRDS(file = "Figure_1/Figure_1d-e/DF_BoxplotsNumberCoefs_Noise.LM22.Figura1.d.RDS")
colnames(LM22_betas_noise)[3] <- "Method"
LM22_betas_noise$Method<-gsub("QUANTISEQ","quanTIseq",LM22_betas_noise$Method) %>% factor(levels=c("ABBAS","ABIS","quanTIseq","CIBERSORT","MIXTURE"))

#Load and format matrix with the output for each algorithm using TIL10 signature of real and estimated number of cell types
TIL10_betas_noise <- readRDS(file = "Figure_1/Figure_1d-e/DF_BoxplotsNumberCoefs_Noise.TIL10.Figura1.d.RDS")
colnames(TIL10_betas_noise)[3] <- "Method"
TIL10_betas_noise$Method<-gsub("QUANTISEQ","quanTIseq",TIL10_betas_noise$Method) %>% factor(levels=c("ABBAS","ABIS","quanTIseq","CIBERSORT","MIXTURE"))

#Plot results for LM22 signature
bx_LM22 <- ggplot(LM22_betas_noise, aes(x=beta, y=est)) +
  geom_boxplot(aes(fill = Method),show.legend=TRUE) +
  labs( x = "True number of coefficients",y = "Estimated number of coefficients", fill = "Method") + 
  scale_y_continuous(name ="Estimated number of coefficients", breaks = c(0:22),labels=as.character(c(0:22))) +
  theme(legend.key = element_blank(),plot.title = element_text(hjust=0.5,size = 16),axis.line = element_line(colour = "black"),
        panel.background = element_blank(),axis.text = element_text(size = 12),axis.title = element_text(size = 14),
        legend.position = "bottom") + 
  ggtitle(label="LM22")
 
#Plot results for TIL10 signature
bx_TIL10 <- ggplot(TIL10_betas_noise, aes(x=beta, y=est)) +
  geom_boxplot(position = position_dodge(width = 0.7), aes(fill = Method),show.legend=TRUE) +
  labs( x = "True number of coefficients",y = "Estimated number of coefficients", fill = "Method") +
  scale_y_continuous(name ="Estimated number of coefficients", breaks = c(0:10),labels=as.character(c(0:10))) +
  theme(plot.title = element_text(hjust=0.5,size = 16),axis.line = element_line(colour = "black"),
        panel.background = element_blank(),axis.text = element_text(size = 12),axis.title = element_text(size = 14),
        legend.position = "bottom") + 
  ggtitle(label="TIL10")

#Save plot in PDF 16x6
ggsave(filename = "Figure_1/Figure_1d-e/Figure_1d-e_celltypes.pdf",plot = grid.arrange(bx_LM22,bx_TIL10,ncol=2),device = "pdf",height = 6,width = 16)

