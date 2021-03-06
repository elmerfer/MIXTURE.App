#Figure_5a and Figure_5h: Immune infiltrate analysis on melanoma biopsies of patients treated with immunotherapy from Van Allen cohort
rm(list=ls())
library(readxl)
library(ggpubr)
library(tibble)
library(gridExtra)
library(grid)
library(dplyr)

#Load MIXTURE output for Van Allen cohort biopsies
load("Figure_5/Figure_5a-h/TEST_MIXTURE_FILE_LM22_vanallen.RData")
mixture_output_vanallen <- mix.test.allen$Subjects$MIXprop

#Load and prepare clinical data
clinical_data_vanallen_1 <- read_excel("Figure_5/Figure_5a-h/clinical_data_vanallen_1.xlsx")
clinical_data_vanallen_2 <- read_excel("Figure_5/Figure_5a-h/clinical_data_vanallen_2.xlsx")
clinical_data_vanallen <- merge(clinical_data_vanallen_1,clinical_data_vanallen_2,by="patient",all.x=TRUE)
#Select features of interest
clinical_data_vanallen<-clinical_data_vanallen[,c(1,4,5,6,17,30,35,59,60,61,63)]
rownames(clinical_data_vanallen)<-clinical_data_vanallen$patient
clinical_data_vanallen<-clinical_data_vanallen[,-c(1,2)]
#Rename features
colnames(clinical_data_vanallen)[c(1,2,3,4,5)]<-c("OS","PFS","group","Estado","RECIST")
clinical_data_vanallen$RECIST[c(14,41)]<-c("PD","SD")

#Calculate Tumor Mutational Burden (as log_nonsynonymous mutation)
for (i in c(1:42)){clinical_data_vanallen$log_nonsyn_mut[i]=log(clinical_data_vanallen$nonsynonymous[i])}
#Split group in "High" or "Low" by the median of TMB
for(i in c(1:42)){if(is.na(clinical_data_vanallen$log_nonsyn_mut[i])==TRUE){clinical_data_vanallen$TMB_group[i]=NA
}else{
  if(clinical_data_vanallen$log_nonsyn_mut[i]<median(na.omit(clinical_data_vanallen$log_nonsyn_mut))){
    clinical_data_vanallen$TMB_group[i]="Low"
  }else{
    clinical_data_vanallen$TMB_group[i]="High"
  }
}
}
#Add Absolute Infiltrate Score for each patient
clinical_data_vanallen<-cbind(clinical_data_vanallen,data.frame(apply(t(mix.test.allen$Subjects$MIXabs[rownames(clinical_data_vanallen),,drop=FALSE]),2,sum)))
colnames(clinical_data_vanallen)[12]<-c("Absolute")
#Classify patients in Responders and Non-responders (pt29 doesnt have RECIST response, it was classified according the "group" feature as Responder)
for (i in c(1:42)){if(clinical_data_vanallen$RECIST[i]=="CR" | clinical_data_vanallen$RECIST[i]=="PR" | clinical_data_vanallen$RECIST[i]=="SD"){clinical_data_vanallen$RECIST_bin[i]="Responder"
}else{
  if(clinical_data_vanallen$RECIST[i]=="PD"){clinical_data_vanallen$RECIST_bin[i]="Non responder"
  }else{clinical_data_vanallen$RECIST_bin[i]="X"}
}
} 
clinical_data_vanallen$RECIST_bin[18]<-c("Responder")

#Merge proportions and clinical datasets
complete_dataset_vanallen <- data.frame()
for (i in c(1:42)){
  for(j in c(1:22)){
    data<-data.frame(rownames = colnames(t(mixture_output_vanallen))[i],Cell_type=rownames(t(mixture_output_vanallen))[j],Proportion=t(mixture_output_vanallen)[j,i])
    complete_dataset_vanallen<-rbind(complete_dataset_vanallen,data)
  }
} 

complete_dataset_vanallen <- complete_dataset_vanallen[order(complete_dataset_vanallen$Cell_type),]
rownames(complete_dataset_vanallen) <- 1:924
complete_dataset_vanallen <- merge(complete_dataset_vanallen,rownames_to_column(clinical_data_vanallen,"rownames"),by="rownames")

#Cell types barplots
#Define y axis digits
scaleFUN <- function(x){sprintf("%.2f", x)}

##Figure_5a: Immune Infiltrate Absolute Score for TMB groups
#Plot
fig5_a1<-ggbarplot(filter(complete_dataset_vanallen,Cell_type=="T cells CD8")[-c(14,41),], x="TMB_group", y="Proportion", color="TMB_group", add = (c("mean_se","jitter"))) + theme(legend.position="none",axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18),axis.title.y = element_text(size = 20)) + scale_color_manual("TMB_group",values=c("palegreen4","yellow3")) + stat_compare_means(method="wilcox.test",comparisons = list(c("High","Low")),bracket.size = 0.5) + xlab("") + ylab("T cells CD8 proportion") + scale_y_continuous(labels=scaleFUN)
fig5_a1[["layers"]][[4]][["aes_params"]]$textsize<-7
fig5_a2<-ggbarplot(filter(complete_dataset_vanallen,Cell_type=="Macrophages M1")[-c(14,41),], x="TMB_group", y="Proportion", color="TMB_group", add = (c("mean_se","jitter"))) + theme(legend.position="none",axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18),axis.title.y = element_text(size = 20))  + scale_color_manual("TMB_group",values=c("palegreen4","yellow3")) + stat_compare_means(method="wilcox.test",comparisons = list(c("High","Low")),label="p.signif") + xlab("") + ylab("Macrophages M1 proportion") + scale_y_continuous(labels=scaleFUN)
fig5_a2[["layers"]][[4]][["aes_params"]]$textsize<-9
fig5_a3<-ggbarplot(filter(complete_dataset_vanallen,Cell_type=="T cells follicular helper")[-c(14,41),], x="TMB_group", y="Proportion", color="TMB_group", add = (c("mean_se","jitter"))) + theme(legend.position="none",axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18),axis.title.y = element_text(size = 20))  + scale_color_manual("TMB_group",values=c("palegreen4","yellow3")) + stat_compare_means(method="wilcox.test",comparisons = list(c("High","Low")),label="p.signif") + xlab("") + ylab("T Cells Follicular Helper proportion") + scale_y_continuous(labels=scaleFUN)
fig5_a3[["layers"]][[4]][["aes_params"]]$textsize<-9

#Save PDF 12x10
ggsave(filename = "Figure_5/Figure_5a-h/Figure_5a_VA.pdf",plot = grid.arrange(ncol=3,nrow=1,fig5_a1,fig5_a2,fig5_a3,top=textGrob("TMB",gp=gpar(fontsize=22))),device = "pdf",width = unit(12,"cm"),height = unit(10,"cm"))

##Figure_5h: Proportion of patients with B cells näive positive infiltrate
#Generate a dataframe with 
table(list(filter(complete_dataset_vanallen,Cell_type=="B\ cells\ naive")$Proportion,filter(complete_dataset_vanallen,Cell_type=="B\ cells\ naive")$RECIST_bin))

#Those with a proportion of 0 are considered "Negative infiltrate, B-", those with proportions higher than 0 are considered "Positive infiltrate,B+"
Bcells_proportion<-data.frame("NR"=c(5,22),"R"=c(8,7))
rownames(Bcells_proportion)<-c("B+","B-")
chisq.test(Bcells_proportion,simulate.p.value = TRUE)

#generate a dataframe with proportions of patients for each group and p-value from the Chi square test
Bcells_proportion<-data.frame("Response"=rep(c("Responder (n=15)","Non responder (n=27)"),2),"B cells naive infiltrate"=c(rep("Positive",2),rep("Negative",2)),"Proportion"=c(Bcells_proportion[1,2]/sum(Bcells_proportion$R),Bcells_proportion[1,1]/sum(Bcells_proportion$NR),Bcells_proportion[2,2]/sum(Bcells_proportion$R),Bcells_proportion[2,1]/sum(Bcells_proportion$NR)))
Bcells_proportion$Response<-factor(Bcells_proportion$Response,levels=c("Responder (n=15)","Non responder (n=27)"))

#Plot and save PDF 9x6
plot_prop_bcells<-ggbarplot(Bcells_proportion,"Response","Proportion", fill="B.cells.naive.infiltrate",palette=c("orange3","royalblue4")) + labs(fill="") + ylab(paste("Proportion of patients with","\nB cells naive infiltrate")) + geom_signif(annotations = c("*"),xmin="Non responder (n=27)",xmax = "Responder (n=15)", y_position = 1.05) + theme(axis.title = element_text(size = 18),legend.text = element_text(size = 16),axis.text = element_text(size=14))
plot_prop_bcells[["layers"]][[2]][["aes_params"]]$textsize<-9
ggsave(filename = "Figure_5/Figure_5a-h/Figure_5h_VA.pdf",plot = plot_prop_bcells,width=unit(6,"cm"),height=unit(9,"cm"))





