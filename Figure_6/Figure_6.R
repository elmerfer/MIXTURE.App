#Figure_6a-d: Analysis of the effect of Anti-CTLA4 treatment prior to Anti-PD1 treatment on a pool of 3 different cohorts
rm(list=ls())
library(tibble)
library(dplyr)
library(ComplexHeatmap)
library(viridis)
library(circlize)
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

##Load Auslander cohort data
load("Figure_5/Figure_5b-c/Auslander_pre_data.Rdata")
#Mixture output
Auslander_pre_proportions <- Auslander_pre_data[,c(8:29)]
rownames(Auslander_pre_proportions) <- paste0("A.",rownames(Auslander_pre_proportions))
#Clinical data
Auslander_pre_clinical_data <- Auslander_pre_data[,-c(8:29)]
rownames(Auslander_pre_clinical_data) <- paste0("A.",rownames(Auslander_pre_clinical_data))

##Merge all datasets
Pool_pre_proportions <- rbind(Riaz_pre_proportions,Liu_pre_proportions,Auslander_pre_proportions)
colnames(Pool_pre_proportions) <- gsub("\\."," ",colnames(Pool_pre_proportions))
Pool_pre_clinical_data <- rbind(Riaz_pre_clinical_data[,c(1,12),drop=FALSE],Liu_pre_clinical_data[,c(46,35),drop=FALSE],Auslander_pre_clinical_data[,c(6,7),drop=FALSE])
data_set_pre_pool <- merge(rownames_to_column(Pool_pre_proportions,"rownames"),rownames_to_column(Pool_pre_clinical_data,"rownames"),by="rownames") %>% column_to_rownames("rownames")
data_set_pre_pool$Response <- factor(data_set_pre_pool$Response,levels = c("Responder","Non responder"))

##Split datasets in two groups: Ipi-naive if biopsies were not treated with Anti-CTLA4 therapy (ipilimumab) and Ipi-progressive if biopsies were treated with Anti_CTLA4 therapy(ipilimumab) and progressed
#Clinical data
Pool_pre_naive_clinical_data <- rownames_to_column(Pool_pre_clinical_data,"rownames") %>% filter(priorCTLA4=="Ipi-naive") %>% column_to_rownames("rownames")
Pool_pre_progressive_clinical_data <- rownames_to_column(Pool_pre_clinical_data,"rownames") %>% filter(priorCTLA4=="Ipi-progressive") %>% column_to_rownames("rownames")

###Ipi-naive group
##Add patients information for each feature
#priorCTLA4
table(Pool_pre_naive_clinical_data$priorCTLA4)
Pool_pre_naive_clinical_data$priorCTLA4 <- gsub("\\Ipi-naive","Ipi-naive (n=102)",Pool_pre_naive_clinical_data$priorCTLA4)
Pool_pre_naive_clinical_data$priorCTLA4 <- factor(Pool_pre_naive_clinical_data$priorCTLA4)
#Add "Study" feature
table(substring(rownames(Pool_pre_naive_clinical_data),1,2))
Pool_pre_naive_clinical_data$Study <- c(rep("Riaz (n=23)",23),rep("Liu (n=73)",73),rep("Auslander (n=6)",6))
Pool_pre_naive_clinical_data$Study <- factor(Pool_pre_naive_clinical_data$Study,levels=c("Auslander (n=6)","Riaz (n=23)","Liu (n=73)"))
#Response
table(Pool_pre_naive_clinical_data$Response)
Pool_pre_naive_clinical_data$Response <- gsub("\\Responder","Responder (n=47)",gsub("\\Non responder","Non responder (n=55)",Pool_pre_naive_clinical_data$Response))
Pool_pre_naive_clinical_data$Response <- factor(Pool_pre_naive_clinical_data$Response,levels = c("Responder (n=47)","Non responder (n=55)"))

###Ipi-progressive
##Add patients information for each feature
#priorCTLA4
table(Pool_pre_progressive_clinical_data$priorCTLA4)
Pool_pre_progressive_clinical_data$priorCTLA4 <- gsub("\\Ipi-progressive","Ipi-progressive (n=59)",Pool_pre_progressive_clinical_data$priorCTLA4)
Pool_pre_progressive_clinical_data$priorCTLA4 <- factor(Pool_pre_progressive_clinical_data$priorCTLA4)
#Add "Study" feature
table(substring(rownames(Pool_pre_progressive_clinical_data),1,2))
Pool_pre_progressive_clinical_data$Study <- c(rep("Riaz (n=26)",26),rep("Liu (n=30)",30),rep("Auslander (n=3)",3))
Pool_pre_progressive_clinical_data$Study <- factor(Pool_pre_progressive_clinical_data$Study,levels=c("Auslander (n=3)","Riaz (n=26)","Liu (n=30)"))
#Response
table(Pool_pre_progressive_clinical_data$Response)
Pool_pre_progressive_clinical_data$Response <- gsub("\\Responder","Responder (n=19)",gsub("\\Non responder","Non responder (n=40)",Pool_pre_progressive_clinical_data$Response))
Pool_pre_progressive_clinical_data$Response <- factor(Pool_pre_progressive_clinical_data$Response,levels = c("Responder (n=19)","Non responder (n=40)"))

###Heatmap
##Generate palette
col_pal<-colorRamp2(seq(-2,2,by=0.02),viridis(length(seq(-2,2,by=0.02))))
##Generate annotations
#Ipi-naive group
pre_pool_naive_top_anno<-HeatmapAnnotation(df=Pool_pre_naive_clinical_data[,c(2,3,1)],which="column",annotation_name_gp = gpar(fontface="bold"),col=list(Response=c("Responder (n=47)"="blue","Non responder (n=55)"="red"),Study=c("Auslander (n=6)"="green4","Riaz (n=23)"="gold1","Liu (n=73)"="tomato4"),priorCTLA4=c("Ipi-naive (n=102)"="navyblue")),gap = unit(1,"mm"),border = TRUE,show_annotation_name = TRUE)
pre_pool_naive_top_anno@anno_list[["priorCTLA4"]]@name<-"Anti CTLA-4 group"
pre_pool_naive_top_anno@anno_list[["priorCTLA4"]]@color_mapping@name<-"Anti CTLA-4 group"
#Ipi-progressive group
pre_pool_progressive_top_anno<-HeatmapAnnotation(df=Pool_pre_progressive_clinical_data[,c(2,3,1)],which="column",annotation_name_gp = gpar(fontface="bold"),col=list(Response=c("Responder (n=19)"="blue","Non responder (n=40)"="red"),Study=c("Auslander (n=3)"="green4","Riaz (n=26)"="gold1","Liu (n=30)"="tomato4"),priorCTLA4=c("Ipi-progressive (n=59)"="darkgreen")),gap = unit(1,"mm"),border = TRUE,show_annotation_name = TRUE)
pre_pool_progressive_top_anno@anno_list[["priorCTLA4"]]@name<-"Anti CTLA-4 group"
pre_pool_progressive_top_anno@anno_list[["priorCTLA4"]]@color_mapping@name<-"Anti CTLA-4 group"

##Plots
#Ipi-naive group (Figure 6a)
HM_pre_pool_naive <- Heatmap(t(scale(Pool_pre_proportions))[,rownames(Pool_pre_naive_clinical_data)],name="Z-score",
                         show_heatmap_legend = TRUE,
                         col=col_pal,border="black",
                         show_column_dend = FALSE,show_row_dend = TRUE,
                         show_column_names = FALSE,show_row_names=TRUE,
                         show_parent_dend_line = FALSE,
                         cluster_rows = TRUE,cluster_columns = TRUE,
                         column_gap = unit(1,"mm"),row_gap = unit(0,"mm"),
                         column_title_gp = gpar(fontsize=12),row_names_gp = gpar(fontsize=10), 
                         top_annotation = pre_pool_naive_top_anno,
                         clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2")

#Save in PDF 13x6
pdf("Figure_6/heatmap_naive.pdf",width=13,height=6)
print(HM_pre_pool_naive)
dev.off()

#Ipi-progressive group (Figure 6b)
HM_pre_pool_progressive <- Heatmap(t(scale(Pool_pre_proportions))[,rownames(Pool_pre_progressive_clinical_data)],name="Z-score",
                               show_heatmap_legend = TRUE,
                               col=col_pal,border="black",
                               show_column_dend = FALSE,show_row_dend = TRUE,
                               show_column_names = FALSE,show_row_names=TRUE,
                               show_parent_dend_line = FALSE,
                               cluster_rows = TRUE,cluster_columns = TRUE,
                               column_gap = unit(1,"mm"),row_gap = unit(0,"mm"),
                               column_title_gp = gpar(fontsize=12),row_names_gp = gpar(fontsize=10), 
                               top_annotation = pre_pool_progressive_top_anno,
                               clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2")

#Save in PDF 13x6
pdf("Figure_6/heatmap_progressive.pdf",width=13,height=6)
print(HM_pre_pool_progressive)
dev.off()

##Cell type proportions statistical analysis: Responders vs Non responders
#Define statistical comparison
response_comp = list(c("Responder","Non responder"))
priorCTLA4_comp = list(c("Ipi-naive","Ipi-progressive"))
#Define y axis digits
scaleFUN <- function(x){sprintf("%.2f", x)}

#Response comparison in groups Ipi-naive and Ipi-treated (Figure 6c)
data_set_pre_pool_aux <- data_set_pre_pool
colnames(data_set_pre_pool_aux) <- gsub(" ",".",colnames(data_set_pre_pool_aux))

BP_CD8_response <- ggbarplot(cbind(data_set_pre_pool_aux[,c(4,23,24)],Cell_type=rep("T cells CD8",161)),x="Response",y="T.cells.CD8",color="Response",facet.by = c("priorCTLA4","Cell_type"),add=(c("mean_se","jitter"))) + theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size = 16),strip.text = element_text(size = 20))+ scale_y_continuous(labels=scaleFUN) +scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method = "wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("T cells CD8 proportion")
BP_CD8_response[["layers"]][[4]][["aes_params"]]$textsize<-8
BP_M2_response <- ggbarplot(cbind(data_set_pre_pool_aux[,c(16,23,24)],Cell_type=rep("Macrophages M2",161)),x="Response",y="Macrophages.M2",color="Response",facet.by = c("priorCTLA4","Cell_type"),add=(c("mean_se","jitter"))) + theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size = 16),strip.text = element_text(size = 20))+ scale_y_continuous(labels=scaleFUN) +scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method = "wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("Macrophages M2 proportion")
BP_M2_response[["layers"]][[4]][["aes_params"]]$textsize<-8
BP_TMR_response <- ggbarplot(cbind(data_set_pre_pool_aux[,c(6,23,24)],Cell_type=rep("T cells CD4 memory resting",161)),x="Response",y="T.cells.CD4.memory.resting",color="Response",facet.by = c("priorCTLA4","Cell_type"),add=(c("mean_se","jitter"))) + theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size = 16),strip.text = element_text(size = 20))+ scale_y_continuous(labels=scaleFUN) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method = "wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("T cells CD4 memory resting proportion")
BP_TMR_response[["layers"]][[4]][["aes_params"]]$textsize<-8
BP_TGD_response <- ggbarplot(cbind(data_set_pre_pool_aux[,c(10,23,24)],Cell_type=rep("T cells gamma delta",161)),x="Response",y="T.cells.gamma.delta",color="Response",facet.by = c("priorCTLA4","Cell_type"),add=(c("mean_se","jitter"))) + theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size = 16),strip.text = element_text(size = 20))+ scale_y_continuous(labels=scaleFUN) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method = "wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("T cells CD4 gamma delta proportion")
BP_TGD_response[["layers"]][[4]][["aes_params"]]$textsize<-8
BP_BN_response <- ggbarplot(cbind(data_set_pre_pool_aux[,c(1,23,24)],Cell_type=rep("B cells naive",161)),x="Response",y="B.cells.naive",color="Response",facet.by = c("priorCTLA4","Cell_type"),add=(c("mean_se","jitter"))) + theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size = 16),strip.text = element_text(size = 20))+ scale_y_continuous(labels=scaleFUN) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method = "wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("B cells naive proportion")
BP_BN_response[["layers"]][[4]][["aes_params"]]$textsize<-8
BP_BM_response <- ggbarplot(cbind(data_set_pre_pool_aux[,c(2,23,24)],Cell_type=rep("B cells memory",161)),x="Response",y="B.cells.memory",color="Response",facet.by = c("priorCTLA4","Cell_type"),add=(c("mean_se","jitter"))) + theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size = 16),strip.text = element_text(size = 20))+ scale_y_continuous(labels=scaleFUN) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method = "wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("B cells naive proportion")
BP_BM_response[["layers"]][[4]][["aes_params"]]$textsize<-8

#Save in PDF 32x20
ggsave(filename = "Figure_6/BP_Ipi_groups_response.pdf",plot=grid.arrange(ncol=6,BP_CD8_response,BP_M2_response,BP_TMR_response,BP_TGD_response,BP_BN_response,BP_BM_response),device = "pdf",width = unit(32,"cm"),height = unit(20,"cm"))

#Responders vs Non responders comparison (Figure 6d)
BP_CD8_priorCTLA4 <- ggbarplot(cbind(data_set_pre_pool_aux[,c(4,23,24)],Cell_type=rep("T cells CD8",161)),x="priorCTLA4",y="T.cells.CD8",color="Response",facet.by = c("Response","Cell_type"),add=(c("mean_se","jitter"))) + theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size = 16),strip.text = element_text(size = 20))+ scale_y_continuous(labels=scaleFUN) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method = "wilcox.test",comparisons = priorCTLA4_comp,label="p.signif") + xlab("") + ylab("T cells CD8 proportion")
BP_CD8_priorCTLA4[["layers"]][[4]][["aes_params"]]$textsize<-8
BP_M1_priorCTLA4 <- ggbarplot(cbind(data_set_pre_pool_aux[,c(15,23,24)],Cell_type=rep("Macrophages M1",161)),x="priorCTLA4",y="Macrophages.M1",color="Response",facet.by = c("Response","Cell_type"),add=(c("mean_se","jitter"))) + theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size = 16),strip.text = element_text(size = 20))+ scale_y_continuous(labels=scaleFUN) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method = "wilcox.test",comparisons = priorCTLA4_comp,label="p.signif") + xlab("") + ylab("Macrophages M1 proportion")
BP_M1_priorCTLA4[["layers"]][[4]][["aes_params"]]$textsize<-8
BP_PC_priorCTLA4 <- ggbarplot(cbind(data_set_pre_pool_aux[,c(3,23,24)],Cell_type=rep("Plasma Cells",161)),x="priorCTLA4",y="Plasma.cells",color="Response",facet.by = c("Response","Cell_type"),add=(c("mean_se","jitter"))) + theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size = 16),strip.text = element_text(size = 20))+ scale_y_continuous(labels=scaleFUN) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method = "wilcox.test",comparisons = priorCTLA4_comp,label="p.signif") + xlab("") + ylab("Plasma Cells proportion")
BP_PC_priorCTLA4[["layers"]][[4]][["aes_params"]]$textsize<-8

#Save in PDF 16x20
ggsave(filename = "Figure_6/BP_Ipi_groups_priorCTLA4.pdf",plot=grid.arrange(ncol=3,BP_CD8_priorCTLA4,BP_M1_priorCTLA4,BP_PC_priorCTLA4),device = "pdf",width = unit(16,"cm"),height = unit(20,"cm"))

