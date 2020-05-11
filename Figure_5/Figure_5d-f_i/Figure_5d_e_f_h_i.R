# Fig 5 d,e,f,i.

library(dplyr)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(limma)

sgm <- read.csv("fpkmvalues_exp.csv") #FPKM generated with count data and pipeline used by Riaz et al. 2017: https://github.com/riazn/bms038_analysis
rownames <- sgm[,1]
sgm <- as.matrix(sgm)
rownames(sgm) <- rownames
sgm <- sgm[,-1]
class(sgm) <- "numeric"
sgm <- avereps[sgm]


# MIXTURE
source("../../Utils/RUN_MIXTURE.R") #Check MIXTURE folder
mixture_riaz_rel <-  mix.test$Subjects$MIXprop
mixture_riaz_rel <- t(mixture_riaz_rel)
mixture_riaz_rel <- as.data.frame(mixture_riaz_rel)

#Clinical data

clinical_data <- read.csv("./clinical_data.csv") # From supplementary information of Riaz et al., 2017. 
patients_groups_riaz <- read.csv("./patients_group.csv") #Table with patiens grouped by time of biopsy
patients.r_riaz <- as.vector(patients_groups_riaz$responders_on)
patients.nr_riaz <- as.vector(patients_groups_riaz$non_responders_on)
patients.nr_riaz <- patients.nr_riaz[patients.nr_riaz!=""]
patients.nr_riaz_notpt32 <- patients.nr_riaz[patients.nr_riaz != "Pt32_On"]

# Filter mixture table for on-treatment biopsies

mixture_matrix_on <- mixture_riaz_rel[,c(patients.r_riaz,patients.nr_riaz_notpt32)]
names(mixture_matrix_on) <- gsub("_.*","",names(mixture_matrix_on))

# Merge and filter clinical data

response_data <- read.table("./SampleTableCorrected.9.19.16.csv", sep=",",header = T,row.names = 1) #Table from https://github.com/riazn/bms038_analysis/tree/master/data
response_data_on <- filter(response_data, PreOn=="On")
clinical_data_on <- merge(response_data_on,clinical_data,by.x="PatientID",by.y="Patient",all.x = T)
clinical_data_on <- filter(clinical_data_on,PatientID %in% names(mixture_matrix_on))

rownames(clinical_data_on) <- clinical_data_on[,1]
clinical_data_on <- clinical_data_on[,c(10,24)]
clinical_data_on$Response.x <- as.character(clinical_data_on$Response.x)
clinical_data_on$response <- clinical_data_on$Response.x
clinical_data_on$response[clinical_data_on$response == "PD"] <- "Non responder"
clinical_data_on$response[clinical_data_on$response == "PRCR" | clinical_data_on$response == "SD"] <- "Responder"
clinical_data_on$Response.x[clinical_data_on$Response.x == "PRCR"] <- "CR/PR"
clinical_data_on$Response.x <- as.factor(clinical_data_on$Response.x)
clinical_data_on$response <- as.factor(clinical_data_on$response)
names(clinical_data_on)[2] <- "Cytolytic score"
names(clinical_data_on)[1] <- "Response"

# Heatmap

mixture_matrix_on_heatmap_scaled <- t(scale(t(mixture_matrix_on)))

col_cyt = colorRamp2(c(0, 600, 600,1700), c("#FFFFFFFF", "#F4C2E3FF", "#E483C8FF", "#DB63BBFF"))

mixture_matrix_on_heatmap_scaled <- mixture_matrix_on_heatmap_scaled[,rownames(clinical_data_on)]
mixture_matrix_on_heatmap_scaled[5,] <- 0
mixture_matrix_on_heatmap_scaled[21,] <- 0

levels(clinical_data_on$response)[levels(clinical_data_on$response)=="Non responder"] <- "Non responder (n=23)"
levels(clinical_data_on$response)[levels(clinical_data_on$response)=="Responder"] <- "Responder (n=31)"

ha = HeatmapAnnotation(Response=clinical_data_on$response,iRecist = clinical_data_on$Response, 
                       `Cytolytic score` = clinical_data_on$`Cytolytic score`,
                       col=list(iRecist=c("CR/PR" = "dodgerblue2", "SD" = "cyan2", "PD"="indianred2"),
                                Response=c(`Responder (n=31)`="blue",`Non responder (n=23)` ="red"),`Cytolytic score`=col_cyt),
                       annotation_name_gp = gpar(fontface="bold"),gap = unit(1,"mm"),na_col = "gray49", border = TRUE,which="column")

color_fig_final <- colorRamp2(seq(-2,2,by=0.02),viridis(length(seq(-2,2,by=0.02))))

plot_heatmap<-Heatmap(mixture_matrix_on_heatmap_scaled,name="Z-score",col=color_fig_final,
                      border="black",show_column_names = FALSE,show_column_dend = T,show_row_names=TRUE, show_row_dend = TRUE,
                      row_names_gp = gpar(fontsize=10), show_parent_dend_line = FALSE,row_gap = unit(0,"mm"),show_heatmap_legend = TRUE,
                      cluster_rows = TRUE,cluster_columns = TRUE,column_gap = unit(1,"mm"),column_title_gp = gpar(fontsize=12),
                      top_annotation = ha,
                      clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2")
draw(plot_heatmap,merge_legend = T) ##Export PDF width=10, height=6


# Paired analysis: pre-on treatment

## Patients paired
patients_groups_riaz_paired <- read.csv("./patients_group_paired.csv") #Table with patiens paired and grouped by time of biopsy
patients.r_riaz_pre_paired <- as.vector(patients_groups_riaz_paired$responders_pre)
patients.r_riaz_on_paired <- as.vector(patients_groups_riaz_paired$responders_on)


## Prpare table
mixture_matrix_pairedSamples_r <- mixture_riaz_rel[,c(patients.r_riaz_pre_paired,patients.r_riaz_on_paired)]
mixture_matrix_pairedSamples_r <- t(mixture_matrix_pairedSamples_r)
mixture_matrix_pairedSamples_r <- as.data.frame(mixture_matrix_pairedSamples_r)
mixture_matrix_pairedSamples_r$Time <- c(rep("Pre",length(patients.r_riaz_pre_paired)),rep("On",length(patients.r_riaz_on_paired)))

## Plots

### M2

plot_m2_r <-  ggpaired(mixture_matrix_pairedSamples_r, x = "Time", y = "`Macrophages M2`", title = "Macrophages M2 - Responders",
                       color = "Time", line.color = "gray", line.size = 1, xlab = F,ylab = "Proportion",
                       palette = "aaas",short.panel.labs = T,combine=F)+theme(text = element_text(size = 12))  +
  stat_compare_means(label = "p.signif", method = "wilcox.test", paired = TRUE, label.x = 1.5,size=5)


### CD4 memory activated

plot_cd4Act_r <-  ggpaired(mixture_matrix_pairedSamples_r, x = "Time", y = "`T cells CD4 memory activated`", title = "T cells CD4 memory activated - Responders",
                           color = "Time", line.color = "gray", line.size = 1, xlab = F,ylab = "Proportion", 
                           palette = "aaas",short.panel.labs = T,combine=F)+theme(text = element_text(size = 12)) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", paired = TRUE, label.x = 1.5,size=5)

### T cells gamma delta

plot_gammaDelt_r <- ggpaired(mixture_matrix_pairedSamples_r, x = "Time", y = "`T cells gamma delta`", title = "T cells gamma delta - Responders",
                             color = "Time", line.color = "gray", line.size = 1, xlab = F,ylab = "Proportion", 
                             palette = "aaas", short.panel.labs = T,combine=F) + theme(text = element_text(size = 12)) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", paired = TRUE, label.x = 1.5,size=5)


## Plots R vs NR on-treatment

data_mixture_riaz_on <- t(mixture_matrix_on)
data_mixture_riaz_on <- merge(data_mixture_riaz_on, clinical_data_on,by.x="row.names",by.y="row.names",all.x=T)
colnames(data_mixture_riaz_on)[2:23] <- gsub(" ", ".",colnames(data_mixture_riaz_on)[2:23])


### CD8 T cells

my_comparisons <- list(c("Non responder", "Responder"))
data_mixture_riaz_on$response <- factor(data_mixture_riaz_on$response,levels = c("Responder","Non responder"))

p <- ggbarplot(data_mixture_riaz_on,x = "response", y = "T.cells.CD8", title = "",
               color = "response", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD8 proportion",
               palette = c("blue","red"),short.panel.labs = T,combine=F, add =c("mean_se","jitter"))+theme(text = element_text(size = 20))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend = "none")


### Macrophages M2

data_mixture_riaz_on$response <- factor(data_mixture_riaz_on$response,levels = c("Responder","Non responder"))

p <- ggbarplot(data_mixture_riaz_on,x = "response", y = "Macrophages.M2", title = "",
               color = "response", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M2 proportion",
               palette = c("blue","red"),short.panel.labs = T,combine=F, add =c("mean_se","jitter"))+theme(text = element_text(size = 20))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend = "none")

### B cells naive

p <- ggbarplot(data_mixture_riaz_on,x = "response", y = "B.cells.naive", title = "",
          color = "response", line.color = "gray", line.size = 1, xlab = F,ylab = "B cells naive proportion",
          palette = c("blue","red"),short.panel.labs = T,combine=F, add =c("mean_se","jitter"))+theme(text = element_text(size = 20))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox", paired = FALSE,size=8)
ggpar(p, legend = "none")


## Chi-squared


### Cell types proportions table on-treatment

sum(filter(data_mixture_riaz_on, response == "Responder") %>% select(B.cells.naive) > 0)
sum(filter(data_mixture_riaz_on, response == "Non responder") %>% select(B.cells.naive) > 0)

sum(filter(data_mixture_riaz_on, response == "Responder") %>% select(B.cells.naive) == 0)
sum(filter(data_mixture_riaz_on, response == "Non responder") %>% select(B.cells.naive) == 0)

length(filter(data_mixture_riaz_on, response == "Responder") %>% select(B.cells.naive) %>% unlist())
length(filter(data_mixture_riaz_on, response == "Non responder") %>% select(B.cells.naive) %>% unlist())

proportions_bnaive_on <-data.frame("Responder"=c(5,26),"Non responder"=c(0,23))
rownames(proportions_bnaive_on)<-c("B naive ++","B naive -")
chisq.test(proportions_bnaive_on,simulate.p.value = T)
proportions_bnaive_on_stack<-data.frame("Response"=rep(c("R (n=31)","NR (n=23)"),2),"B cells naive infiltrate"=c(rep("Positive",2),rep("Negative",2)),"Proportion"=c(0.1612903,0,0.8387097,1))



### Stacked plot - PDF 9x6

proportions_bnaive_on_stack$Response <- c("Responders (n=31)", "Non responder (n=23)", "Responders (n=31)", "Non responder (n=23)")

proportions_bnaive_on_stack$Response <- factor(proportions_bnaive_on_stack$Response,levels = c("Responders (n=31)","Non responder (n=23)"))

plot_prop_bcells<-ggbarplot(proportions_bnaive_on_stack,"Response","Proportion", fill="B.cells.naive.infiltrate",palette=c("orange3","royalblue4")) + labs(fill="") + ylab(paste("Proportion of patients with","\nB cells naive infiltrate")) + geom_signif(annotations = c("0.06247"), xmin = "Responders (n=31)",xmax = "Non responder (n=23)", y_position = 1.05) + theme(axis.title = element_text(size = 18),legend.text = element_text(size = 16),axis.text = element_text(size=14))
plot_prop_bcells[["layers"]][[2]][["aes_params"]]$textsize<-5
plot_prop_bcells


