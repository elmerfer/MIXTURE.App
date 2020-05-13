#Figure_5b-c: Immune infiltrate analysis of a pool of biopsies from four cohorts treated with anti-PD1 therapy
rm(list=ls())
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(dplyr)
library(tibble)
library(ggpubr)
library(gridExtra)

##Load Riaz cohort data
load("Figure_5/Figure_5b-c/Riaz_pre_data.Rdata")
#Mixture output
Riaz_pre_proportions <- Riaz_pre_data[,c(1:22)]
rownames(Riaz_pre_proportions) <- paste0("R.",rownames(Riaz_pre_proportions))
#Clinical data
Riaz_pre_clinical_data <- Riaz_pre_data[,-c(1:22)]
rownames(Riaz_pre_clinical_data) <- paste0("R.",rownames(Riaz_pre_clinical_data))

##Load Liu cohort data
load("Figure_5/Figure_5b-c/Liu_pre_data.Rdata")
#Mixture output
Liu_pre_proportions <- Liu_pre_data[,c(46:67)]
rownames(Liu_pre_proportions) <- paste0("L.",rownames(Liu_pre_proportions))
#Clinical data
Liu_pre_clinical_data <- Liu_pre_data[,-c(46:67)]
rownames(Liu_pre_clinical_data) <- paste0("L.",rownames(Liu_pre_clinical_data))
#Classify patients response according to RECIST
for (i in c(1:nrow(Liu_pre_clinical_data))){
  if(Liu_pre_clinical_data$RECIST[i]=="PD"){Liu_pre_clinical_data$Response[i]="Non responder"
  }else{
    Liu_pre_clinical_data$Response[i]="Responder"
  }
}

##Load Hugo cohort data
load("Figure_5/Figure_5b-c/Hugo_pre_data.Rdata")
#Mixture output
Hugo_pre_proportions <- Hugo_pre_data[,c(14:35)]
rownames(Hugo_pre_proportions) <- paste0("H.",rownames(Hugo_pre_proportions))
#Clinical data
Hugo_pre_clinical_data <- Hugo_pre_data[,-c(14:35)]
rownames(Hugo_pre_clinical_data) <- paste0("H.",rownames(Hugo_pre_clinical_data))
#Classify patients response according to RECIST
for (i in c(1:nrow(Hugo_pre_clinical_data))){if (Hugo_pre_clinical_data$irRECIST[i]=="PD"){Hugo_pre_clinical_data$Response[i]="Non responder"
}else{
  Hugo_pre_clinical_data$Response[i]="Responder"}
}

##Load Auslander cohort data
load("Figure_5/Figure_5b-c/Auslander_pre_data.Rdata")
#Mixture output
Auslander_pre_proportions <- Auslander_pre_data[,c(8:29)]
rownames(Auslander_pre_proportions) <- paste0("A.",rownames(Auslander_pre_proportions))
#Clinical data
Auslander_pre_clinical_data <- Auslander_pre_data[,-c(8:29)]
rownames(Auslander_pre_clinical_data) <- paste0("A.",rownames(Auslander_pre_clinical_data))

##Merge all datasets
Pool_pre_proportions <- rbind(Riaz_pre_proportions,Liu_pre_proportions,Hugo_pre_proportions,Auslander_pre_proportions)
colnames(Pool_pre_proportions) <- gsub("\\."," ",colnames(Pool_pre_proportions))
Pool_pre_clinical_data <- rbind(Riaz_pre_clinical_data[,1,drop=FALSE],Liu_pre_clinical_data[,46,drop=FALSE],Hugo_pre_clinical_data[,14,drop=FALSE],Auslander_pre_clinical_data[,6,drop=FALSE])
data_set_pre_pool <- merge(rownames_to_column(Pool_pre_proportions,"rownames"),rownames_to_column(Pool_pre_clinical_data,"rownames"),by="rownames") %>% column_to_rownames("rownames")

#Show number of patients from each cohort
table(substring(rownames(data_set_pre_pool),1,1))
#Add the "Study" feature for each patient
data_set_pre_pool <- cbind(data_set_pre_pool,Study=c(rep("Auslander (n=9)",9),rep("Hugo (n=27)",27),rep("Liu (n=103)",103),rep("Riaz (n=49)",49)))
data_set_pre_pool$Study <- factor(data_set_pre_pool$Study,levels=c("Auslander (n=9)","Hugo (n=27)","Riaz (n=49)","Liu (n=103)"))
#Add number of Responder and Non responder patients
data_set_pre_pool$Response <- gsub("\\Responder","Responder (n=81)",gsub("\\Non responder","Non responder (n=107)",data_set_pre_pool$Response))
data_set_pre_pool$Response <- factor(data_set_pre_pool$Response,levels=c("Responder (n=81)","Non responder (n=107)"))

##Heatmap
#Generate annotation
pre_pool_top_anno <- HeatmapAnnotation(df=data_set_pre_pool[,c(24,23),drop=FALSE],annotation_name_gp = gpar(fontface="bold"),gap = unit(1,"mm"),border = TRUE,which="column",col=list(Response=c("Responder (n=81)"="blue","Non responder (n=107)"="red"),Study=c("Auslander (n=9)"="green4","Hugo (n=27)"="navyblue","Riaz (n=49)"="gold1","Liu (n=103)"="tomato4")))
#Generate palette
col_pal<-colorRamp2(seq(-2,2,by=0.02),viridis(length(seq(-2,2,by=0.02))))
#Plot
plot_pool_pre <- Heatmap(t(scale(data_set_pre_pool[,c(1:22)])),name="Z-score",show_heatmap_legend = TRUE,
                   col=col_pal,border="black",
                   show_column_dend = FALSE,show_row_dend = TRUE,
                   show_column_names = FALSE,show_row_names=TRUE,
                   show_parent_dend_line = FALSE,
                   cluster_rows = TRUE,cluster_columns = TRUE,
                   column_gap = unit(1,"mm"),row_gap = unit(0,"mm"),
                   column_title_gp = gpar(fontsize=12),row_names_gp = gpar(fontsize=10), 
                   top_annotation = pre_pool_top_anno,column_split = data_set_pre_pool$Response,
                   clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2")

pdf("Figure_5/Figure_5b-c/Pool_pre_heatmap.pdf",width=13,height=6)
print(plot_pool_pre)
dev.off()

##Cell type proportions statistical analysis: Responders vs Non responders
#Define statistical comparison
response_comp=list(c("Responder","Non responder"))
#Define y axis digits
scaleFUN <- function(x){sprintf("%.2f", x)}
#Modify response levels to erase number of patients from the plot
data_set_pre_pool_aux <- data_set_pre_pool
data_set_pre_pool_aux$Response <- as.character(data_set_pre_pool_aux$Response)
data_set_pre_pool_aux$Response[which(data_set_pre_pool_aux$Response=="Responder (n=81)")] <- "Responder"
data_set_pre_pool_aux$Response[which(data_set_pre_pool_aux$Response=="Non responder (n=107)")] <- "Non responder"
data_set_pre_pool_aux$Response <- factor(data_set_pre_pool_aux$Response,levels = c("Responder","Non responder"))

colnames(data_set_pre_pool_aux) <- gsub(" ",".",colnames(data_set_pre_pool_aux))
#Plot
BP_pool_BN <- ggbarplot(data_set_pre_pool_aux, x="Response", y="B.cells.naive", color="Response", add = (c("mean_se","jitter"))) + theme(legend.position="none",axis.text = element_text(size = 22),axis.title.y = element_text(size = 22)) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method="wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("B cells Naive proportion") + scale_y_continuous(labels=scaleFUN)
BP_pool_BN[["layers"]][[4]][["aes_params"]]$textsize<-10
BP_pool_BM <- ggbarplot(data_set_pre_pool_aux, x="Response", y="B.cells.memory", color="Response", add = (c("mean_se","jitter"))) + theme(legend.position="none",axis.text = element_text(size = 22),axis.title.y = element_text(size = 22)) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method="wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("B cells Memory proportion") + scale_y_continuous(labels=scaleFUN)
BP_pool_BN[["layers"]][[4]][["aes_params"]]$textsize<-10
BP_pool_CD8 <- ggbarplot(data_set_pre_pool_aux, x="Response", y="T.cells.CD8", color="Response", add = (c("mean_se","jitter"))) + theme(legend.position="none",axis.text = element_text(size = 22),axis.title.y = element_text(size = 22)) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method="wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("T cells CD8 proportion") + scale_y_continuous(labels=scaleFUN)
BP_pool_CD8[["layers"]][[4]][["aes_params"]]$textsize<-10
BP_pool_TMR <- ggbarplot(data_set_pre_pool_aux, x="Response", y="T.cells.CD4.memory.resting", color="Response", add = (c("mean_se","jitter"))) + theme(legend.position="none",axis.text = element_text(size = 22),axis.title.y = element_text(size = 22)) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method="wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("T cells CD4 Memory resting proportion") + scale_y_continuous(labels=scaleFUN)
BP_pool_TMR[["layers"]][[4]][["aes_params"]]$textsize<-10
BP_pool_TGD <- ggbarplot(data_set_pre_pool_aux, x="Response", y="T.cells.gamma.delta", color="Response", add = (c("mean_se","jitter"))) + theme(legend.position="none",axis.text = element_text(size = 22),axis.title.y = element_text(size = 22)) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method="wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("T cells Gamma Delta proportion") + scale_y_continuous(labels=scaleFUN)
BP_pool_TGD[["layers"]][[4]][["aes_params"]]$textsize<-10
BP_pool_M2 <- ggbarplot(data_set_pre_pool_aux, x="Response", y="Macrophages.M2", color="Response", add = (c("mean_se","jitter"))) + theme(legend.position="none",axis.text = element_text(size = 22),axis.title.y = element_text(size = 22)) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method="wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("Macrophages M2 proportion") + scale_y_continuous(labels=scaleFUN)
BP_pool_M2[["layers"]][[4]][["aes_params"]]$textsize<-10

#Save all plots in PDF 32x12
ggsave(filename = "Figure_5/Figure_5b-c/Boxplots_pool_prePD1.pdf",plot = grid.arrange(BP_pool_CD8,BP_pool_M2,BP_pool_TMR,BP_pool_TGD,BP_pool_BN,BP_pool_BM,ncol=6),device = "pdf",width = unit(32,"cm"),height = unit(12,"cm"))



