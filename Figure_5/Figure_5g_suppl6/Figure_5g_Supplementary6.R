#Figure_5g: Immune Absolute Score in Liu cohort and Supplementary_Figure_6.
rm(list=ls())
library(readxl)
library(ggpubr)
library(gridExtra)

##Load mixture output for Liu cohort
mixture_output_absolute_Liu <- read.csv("Figure_5/Figure_5g_suppl6/Immune_Absolute_Score_MIXTURE_Liu.csv",sep=";",row.names = "X")
#Exclude Patient107 who had low gene counts and therefore MIXTURE wasn't able to estimate it's immune infiltrate
mixture_output_absolute_Liu <- mixture_output_absolute_Liu[-c(grep("107",rownames(mixture_output_absolute_Liu))),]

###Load clinical data for Liu cohort
clinical_data_Liu <- read.csv("Figure_5/Figure_5g_suppl6/clinicaldata_Liu.csv",sep=";",row.names="patient")

##Exclude Patients that doesn't have transcriptomic data + Patient107
clinical_data_Liu <- clinical_data_Liu[rownames(mixture_output_absolute_Liu),]

##Filter features of interest
clinical_data_Liu <- clinical_data_Liu[,c(2,5,11)]

##Calculate Absolute Score
clinical_data_Liu <- cbind(clinical_data_Liu,Absolute=apply(mixture_output_absolute_Liu,1,sum))

##Work with Heterogeneity
#Define Hig and Low Heteroegeneity groups
for (i in c(1:nrow(clinical_data_Liu))){
  if(clinical_data_Liu$heterogeneity[i]<median(clinical_data_Liu$heterogeneity)){
    clinical_data_Liu$ITH_group[i]="Low"
  }else{
    clinical_data_Liu$ITH_group[i]="High"}
}
clinical_data_Liu$ITH_group <- factor(clinical_data_Liu$ITH_group)

##Work with number of mutations
#Calculate TMB as the log of nonsynonymous mutation
for (i in c(1:nrow(clinical_data_Liu))){clinical_data_Liu$Log_nonsyn_mut[i]=log(clinical_data_Liu$nonsyn_muts[i])}
#Define High and Low TMB groups
for (i in c(1:nrow(clinical_data_Liu))){
  if(clinical_data_Liu$Log_nonsyn_mut[i]<median(clinical_data_Liu$Log_nonsyn_mut)){
    clinical_data_Liu$TMB_group[i]="Low"
  }else{
      clinical_data_Liu$TMB_group[i]="High"}
}
clinical_data_Liu$TMB_group <- factor(clinical_data_Liu$TMB_group)

##Work with response
for (i in c(1:nrow(clinical_data_Liu))){
  if(clinical_data_Liu$BR[i]=="PD"){clinical_data_Liu$Response[i]="Non responder"
  }else{
    clinical_data_Liu$Response[i]="Responder"}
}
clinical_data_Liu$Response <- factor(clinical_data_Liu$Response,levels=c("Responder","Non responder"))

##Work with TMB and heterogeneity groups
for (i in c(1:nrow(clinical_data_Liu))){
  if(clinical_data_Liu$ITH_group[i]=="High" & clinical_data_Liu$TMB_group[i]=="High"){
    clinical_data_Liu$ITH_TMB_group[i]="High_ITH + High_TMB"
  }else{
    if(clinical_data_Liu$ITH_group[i]=="Low" & clinical_data_Liu$TMB_group[i]=="Low"){
      clinical_data_Liu$ITH_TMB_group[i]="Low_ITH + Low_TMB"
    }else{
      if(clinical_data_Liu$ITH_group[i]=="High" & clinical_data_Liu$TMB_group[i]=="Low"){
        clinical_data_Liu$ITH_TMB_group[i]="High_ITH + Low_TMB"
      }else{
        clinical_data_Liu$ITH_TMB_group[i]="Low_ITH + High_TMB"
      }
    }
  }
}
clinical_data_Liu$ITH_TMB_group <- factor(clinical_data_Liu$ITH_TMB_group,levels=c("High_ITH + High_TMB","Low_ITH + Low_TMB","High_ITH + Low_TMB","Low_ITH + High_TMB"))

##Reorder features
clinical_data_Liu <- clinical_data_Liu[,c(3,8,2,7,1,5,6,4,9)]

##Cell type proportions statistical analysis: Responders vs Non responders
#Define statistical comparison
response_comp=list(c("Responder","Non responder"))
#Define y axis digits
scaleFUN <- function(x){sprintf("%.2f", x)}

#Plot Absolute score for ITH and TMB groups              
plot_ITHTMB_Liu <- ggbarplot(clinical_data_Liu, x="Response", y="Absolute", color="Response",facet.by = c("","ITH_TMB_group"), add = (c("mean_se","jitter"))) + theme(legend.position="none") + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method="wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("Infiltrate Absolute Score") + scale_y_continuous(labels=scaleFUN) + theme(axis.title = element_text(size = 16),strip.text = element_text(size = 18),axis.text = element_text(size = 14))
plot_ITHTMB_Liu[["layers"]][[4]][["aes_params"]]$textsize<-7

#Save all plots in PDF 12x7
ggsave(filename = "Figure_5/Figure_5g_suppl6/Barplot_ITH_TMB.pdf",plot = plot_ITHTMB_Liu,device = "pdf",width = unit(12,"cm"),height = unit(7,"cm"))

#Supplementary Figure 6
supl5b<-ggbarplot(clinical_data_Liu, x="Response", y="Absolute", color="Response",facet.by = "ITH_group", add = (c("mean_se","jitter"))) + theme(legend.position="none",plot.title = element_text(hjust = 0.5) ) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method="wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("Infiltrate Absolute Score")  + scale_y_continuous(labels=scaleFUN) + ggtitle(label="ITH",)
supl5c<-ggbarplot(clinical_data_Liu, x="Response", y="Absolute", color="Response",facet.by = "TMB_group", add = (c("mean_se","jitter"))) + theme(legend.position="none",plot.title = element_text(hjust = 0.5)) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method="wilcox.test",comparisons = response_comp,label="p.signif") + xlab("") + ylab("Infiltrate Absolute Score") + scale_y_continuous(labels=scaleFUN) + ggtitle(label="TMB")

clinical_data_Liu$Response<-gsub("^Responder$","Responder (n=64)",gsub("^Non responder$","Non responder (n=56)",clinical_data_Liu$Response))
clinical_data_Liu$Response<-factor(clinical_data_Liu$Response,levels = c("Responder (n=64)","Non responder (n=56)"))
supl5a<-ggbarplot(clinical_data_Liu, x="Response", y="Absolute", color="Response", add = (c("mean_se","jitter"))) + theme(legend.position="none",) + scale_color_manual("Response",values=c("blue","red")) + stat_compare_means(method="wilcox.test",comparisons = list(c("Responder (n=64)","Non responder (n=56)")),label="p.signif") + xlab("") + ylab("Infiltrate Absolute Score") + scale_y_continuous(labels=scaleFUN)

#Save all plots in PDF 21x5
ggsave(filename = "Figure_5/Figure_5g_suppl6/supp6_Absolute.pdf",plot = grid.arrange(ncol=3,nrow=1,supl5a,supl5b,supl5c),device = "pdf",width = unit(21,"cm"),height = unit(5,"cm"))

