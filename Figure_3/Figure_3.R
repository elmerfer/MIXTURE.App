#Figure_3: Immune infiltrate analysis on TCGA-BRCA dataset
rm(list=ls())
library(ggpubr)
library(gridExtra)
library(RGraphics)
library(ggcorrplot)
library(ComplexHeatmap)

##Cell type proportions analysis for 3 algorithms (quanTIseq, CIBERSORT and MIXTURE) using TIL10 or LM22 signature matrix (Figure_3a)
#Load mixture output for TCGA_BRCA dataset
load(file = "Figure_3//TCGA-BRCA_output_MIXTURE.RData")

#Plot LM22 results
LM22_BRCA <- ggplot(mixture_output_BRCA_LM22, aes(x=CT, y=betahat, fill = Method)) + ylab(label="") + theme_bw() + scale_color_manual("Method",values=c("quanTIseq"="seagreen4","CIBERSORT"="steelblue2","MIXTURE"="mediumorchid2"),aesthetics = "fill") + theme(legend.position = "none",axis.text.x = element_blank(),strip.background.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.text.x = element_text(size=20),axis.text.y = element_text(size = 16)) + geom_boxplot() + ylim(0,1) + facet_grid(~Method) + xlab(label="") + ylab(label="LM22")

#Plot TIL10 results
TIL10_BRCA<-ggplot(mixture_output_BRCA_til10, aes(x=CT, y=betahat, fill = Method)) + ylab(label="") + scale_color_manual("Method",values=c("quanTIseq"="seagreen4","CIBERSORT"="steelblue2","MIXTURE"="mediumorchid2"),aesthetics = "fill") + theme_bw() + geom_boxplot() + theme(axis.title.x = element_text(size = 20),legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1,size = 16),strip.text.x = element_text(size = 0),strip.background.x = element_blank(),axis.text.y = element_text(size = 16)) + facet_grid(~Method) + xlab(label="") + ylab(label = "TIL10")  

#Save plot in PDF 20x12
ggsave(filename = "Figure_3/BP_BRCA.pdf",plot=grid.arrange(LM22_BRCA,TIL10_BRCA, nrow = 2,left=textGrob(label="Estimated cell type proportion",hjust = 0.5,rot = 90,gp = gpar(fontsize=20))),device = "pdf",width=unit(20,"cm"),height = unit(12,"cm"))

##Correlation between cell types proportions estimated by MIXTURE using LM22 or TIL10 signatures (Figure_3b)
LM22_TIL10_cor_mat <- readRDS("Figure_3/LM22_TIL10_correlation_BRCA.RDS")

#Plot
LM22_TIL10_corplot <- ggcorrplot(t(LM22_TIL10_cor_mat), method = "circle")

#Save plot in PDF 3x6
ggsave(filename = "Figure_3/Corplot_LM22_TIL10_BRCA.pdf",plot=LM22_TIL10_corplot,device = "pdf",width=unit(3,"cm"),height=unit(6,"cm"))

##Cox hazar analysis for data using LM22 matrix (Figure_3c)
Method.test <- "MIXTURE"
Signature <- "LM22"

#Load output from Cox Analysis
load(file = "Figure_3/DF_forplot_CoxAnalisisMIXTURE_FPKM_LM22.RData")

#Format matrix
cox_MIXTURE_LM22$Cox <- as.character(cox_MIXTURE_LM22$Cox)
cox_MIXTURE_LM22$Cox[which(is.na(cox_MIXTURE_LM22$Cox))] <- 0
cox_MIXTURE_LM22_cor <- data.frame(row.names = unique(cox_MIXTURE_LM22$ct),All=cox_MIXTURE_LM22$Cox[1:22],LumB=cox_MIXTURE_LM22$Cox[23:44],
                              Basal=cox_MIXTURE_LM22$Cox[45:66],LumA=cox_MIXTURE_LM22$Cox[67:88],Her2=cox_MIXTURE_LM22$Cox[89:110],
                              ER=cox_MIXTURE_LM22$Cox[111:132],ER=cox_MIXTURE_LM22$Cox[133:154],stringsAsFactors = FALSE)
cox_MIXTURE_LM22_cor <- data.matrix(cox_MIXTURE_LM22_cor)
colnames(cox_MIXTURE_LM22_cor)[c(6,7)] <- c("ER+","ER-")
cox_MIXTURE_LM22_cor <- cox_MIXTURE_LM22_cor[seq(from=22,to=1),]

#Plot
cox_BRCA_plot <- ggcorrplot(t(cox_MIXTURE_LM22_cor), method = "circle",title = "MIXTURE",colors = c("green","white","red"),show.legend = FALSE) + 
                 theme(plot.title = element_text(hjust = 0.5,size = 18))
#Generate legend
legend <- legendGrob(labels = c("Prognosis","Good","Poor"),pch = 20,gp = gpar(col=c("white","green","red")))

cox_BRCA_plot_final <- grid.arrange(ncol=1,nrow=1,plot=cox_BRCA_plot,left=textGrob(label="LM22",vjust = 0.6,hjust = 0.2,rot=90,gp = gpar(fontsize=18)),
                                    bottom=textGrob(label="Molecular Subtype",hjust=0.02,gp = gpar(fontsize=14)),right=legend)

#Save plot in PDF 10x8
ggsave(filename = "Figure_3/Cox_plot_BRCA.pdf",plot = cox_BRCA_plot_final,device = "pdf",width=unit(14,"cm"),height=unit(10,"cm"))







