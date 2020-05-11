#Figure 4

#Download data from TCGA

library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)

##SKCM
query_sckm_fpkm <- GDCquery(project = "TCGA-SKCM",
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "HTSeq - FPKM")
GDCdownload(query_sckm_fpkm)
data_skcm_fpkm <- GDCprepare(query_sckm_fpkm)

##HNSC
query_HNSC_fpkm <- GDCquery(project = "TCGA-HNSC",
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "HTSeq - FPKM")
GDCdownload(query_HNSC_fpkm)
data_HNSC_fpkm <- GDCprepare(query_HNSC_fpkm)

##LUAD
query_luad_fpkm <- GDCquery(project = "TCGA-LUAD",
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "HTSeq - FPKM")
GDCdownload(query_luad_fpkm)
data_luad_fpkm <- GDCprepare(query_luad_fpkm)

##LUAD counts
query_luad_counts <- GDCquery(project = c("TCGA-LUAD"), 
                              data.category = "Transcriptome Profiling", 
                              data.type = "Gene Expression Quantification", 
                              workflow.type = "HTSeq - Counts")

GDCdownload(query_luad_counts)
expdat_luad_counts <- GDCprepare(query = query_luad_counts)


##COAD
query_coad_fpkm <- GDCquery(project = "TCGA-COAD",
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "HTSeq - FPKM")
GDCdownload(query_coad_fpkm)
data_coad_fpkm <- GDCprepare(query_coad_fpkm)

##COAD counts
query_coad_counts <- GDCquery(project = "TCGA-COAD",
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "HTSeq - Counts")
GDCdownload(query_coad_counts)
data_coad_counts <- GDCprepare(query_coad_counts)


#Prepare data
##summerize experiment: extract rowData, colData, assay

##MELANOMA-HTSEQ-FPKM
gene_data_melanoma <- rowData(data_skcm_fpkm)
assay_data_melanoma <- assay(data_skcm_fpkm)
colData_melanoma <- colData(data_skcm_fpkm)
rownames(assay_data_melanoma) <- gene_data_melanoma$external_gene_name
assay_data_melanoma_df <- data.frame(assay_data_melanoma, check.names = F)
colnames(assay_data_melanoma_df) <- colData_melanoma$bcr_patient_barcode

##HNSC-HTSEQ-FPKM
gene_data_HNSC <- rowData(data_HNSC_fpkm)
assay_data_HNSC <- assay(data_HNSC_fpkm)
colData_hnsc <- colData(data_HNSC_fpkm)
rownames(assay_data_HNSC) <- gene_data_HNSC$external_gene_name
assay_data_HNSC_df <- data.frame(assay_data_HNSC, check.names = F)
colnames(assay_data_HNSC_df) <- colData_hnsc$bcr_patient_barcode

##LUAD-HTSEQ-FPKM
gene_data_luad <- rowData(data_luad_fpkm)
assay_data_luad <- assay(data_luad_fpkm)
colData_luad <- colData(data_luad_fpkm)
rownames(assay_data_luad) <- gene_data_luad$external_gene_name
assay_data_luad_df <- data.frame(assay_data_luad, check.names = F)
colnames(assay_data_luad_df) <- colData_luad$bcr_patient_barcode

##LUAD-HTSEQ-Counts
gene_data_luad_counts<- rowData(expdat_luad_counts)
assay_data_luad_counts <- assay(expdat_luad_counts)
colData_luad_counts <- colData(expdat_luad_counts)
rownames(assay_data_luad_counts) <- gene_data_luad_counts$external_gene_name
assay_data_luad_counts_df <- data.frame(assay_data_luad_counts, check.names = F)
colnames(assay_data_luad_counts_df) <- colData_luad_counts$bcr_patient_barcode

##COAD-HTSEQ-FPKM
gene_data_COAD <- rowData(data_coad_fpkm)
assay_data_COAD <- assay(data_coad_fpkm)
colData_COAD <- colData(data_coad_fpkm)
rownames(assay_data_COAD) <- gene_data_COAD$external_gene_name
assay_data_COAD_df <- data.frame(assay_data_COAD, check.names = F)
colnames(assay_data_COAD_df) <- colData_COAD$bcr_patient_barcode

##COAD-HTSEQ-Counts
gene_data_COAD_counts <- rowData(data_coad_counts)
assay_data_COAD_counts <- assay(data_coad_counts)
colData_COAD_counts <- colData(data_coad_counts)
rownames(assay_data_COAD_counts) <- gene_data_COAD_counts$external_gene_name
assay_data_COAD_counts_df <- data.frame(assay_data_COAD_counts, check.names = F)
colnames(assay_data_COAD_counts_df) <- colData_COAD_counts$bcr_patient_barcode


#TCGA-LUAD analysis

##Filter normal samples and replicates

tumor_samples_luad <- colnames(assay_data_luad_df)[substr(colnames(assay_data_luad_df),14,15) !=11]
tumor_samples_luad_replicates <- colnames(assay_data_luad_df)[duplicated(substr(colnames(assay_data_luad_df),1,16))]

assay_data_luad_df <- assay_data_luad_df[,tumor_samples_luad]
assay_data_luad_df <- assay_data_luad_df[,!colnames(assay_data_luad_df) %in% tumor_samples_luad_replicates]
assay_data_luad_counts_df <- assay_data_luad_counts_df[, tumor_samples_luad]
assay_data_luad_counts_df <- assay_data_luad_counts_df[,!colnames(assay_data_luad_counts_df) %in% tumor_samples_luad_replicates]


##MIXTURE
sgm <- assay_data_luad_df
source("../../Utils/RUN_MIXTURE.R") #Check RUN_MIXTURE.R location - This run MIXTURE
#mix.test object is generated after run
results_mixture_luad_fpkm <-  mix.test$Subjects$MIXprop #relative proportions
results_mixture_luad_fpkm <- t(results_mixture_luad_fpkm)
results_mixture_luad_fpkm <- as.data.frame(results_mixture_luad_fpkm)


##PDL1 expression

luad_expression_tmm <- DGEList(counts=assay_data_luad_counts_df)
luad_expression_tmm <- calcNormFactors(luad_expression_tmm, method = "TMM")
tmm_luad <- cpm(luad_expression_tmm)
tmm_luad_pdl1 <- tmm_luad["CD274",,drop=F]
tmm_luad_pdl1 <- t(tmm_luad_pdl1)
tmm_luad_pdl1 <- as.data.frame(tmm_luad_pdl1)
tmm_luad_pdl1$Tumor_Sample_Barcode <- rownames(tmm_luad_pdl1)
rownames(tmm_luad_pdl1) <- NULL
tmm_luad_pdl1 <- tmm_luad_pdl1[,c(2,1)]
tmm_luad_pdl1 <- mutate(tmm_luad_pdl1, short_sample_code = substr(tmm_luad_pdl1$Tumor_Sample_Barcode,1,15))

##LUAD - TMB from Cell of Origin proyect (https://gdc.cancer.gov/about-data/publications/PanCan-CellOfOrigin)

tmb_cell_of_origin <- read.table("./mutation-load_updated.txt",header = T,check.names = F,sep = "\t")
tmb_luad <- filter(tmb_cell_of_origin,Tumor_Sample_ID %in% substr(tumor_samples_luad,1,15))
tmb_luad <- mutate(tmb_luad, TMB = tmb_luad$`Silent per Mb` + tmb_luad$`Non-silent per Mb`)
tmb_luad <- tmb_luad[,-c(1,2)]

luad_tmb_quantiles <- quantile(tmb_luad$TMB, probs = c(0.33, 0.66))

luad_tmb_groups <- ifelse(tmb_luad$TMB < luad_tmb_quantiles[1], "Low",
                          ifelse(tmb_luad$TMB >= luad_tmb_quantiles[1] & tmb_luad$TMB < luad_tmb_quantiles[2] , "Intermediate",
                                 ifelse(tmb_luad$TMB >= luad_tmb_quantiles[2], "High", NA)))

tmb_luad$TMB_Group <- as.vector(luad_tmb_groups)
luad_features$TMB_Group <- factor(luad_features$TMB_Group, levels = c("Low","Intermediate","High"))

luad_features <- merge(tmm_luad_pdl1,tmb_luad,by.x="short_sample_code",by.y="Tumor_Sample_ID",all.x=T)

#Mutations

library(TCGAmutations)
library(maftools)

tcga_load("LUAD")

patients_luad_egfr_mut <- filter(tcga_luad_mc3@data,Hugo_Symbol == "EGFR") %>% select(Tumor_Sample_Barcode_full) %>% unique()
patients_luad_l858r_mut <- filter(tcga_luad_mc3@data, HGVSp_Short == "p.L858R" & Hugo_Symbol == "EGFR") %>% select(Tumor_Sample_Barcode_full)
patients_luad_egfr_del19 <- filter(tcga_luad_mc3@data, Variant_Type == "DEL" & Hugo_Symbol == "EGFR")
patients_luad_egfr_del19 <- mutate(patients_luad_egfr_del19, HGVSp_positionInit = substr(HGVSp_Short,4,6))
patients_luad_egfr_del19$HGVSp_positionInit <- as.numeric(patients_luad_egfr_del19$HGVSp_positionInit)
patients_luad_egfr_del19 <- filter(patients_luad_egfr_del19, HGVSp_positionInit > 728 & HGVSp_positionInit < 762) #Exon 19
patients_luad_egfr_del19 <- select(patients_luad_egfr_del19, Tumor_Sample_Barcode_full) %>% unique()
patients_luad_tp53_mut <- filter(tcga_luad_mc3@data,Hugo_Symbol == "TP53") %>% select(Tumor_Sample_Barcode_full) %>% unique()
patients_luad_stk11_mut <- filter(tcga_luad_mc3@data,Hugo_Symbol == "STK11") %>% select(Tumor_Sample_Barcode_full) %>% unique()

luad_egfr_state <- ifelse(luad_features$Tumor_Sample_Barcode %in% substr(patients_luad_egfr_mut$Tumor_Sample_Barcode_full,1,16), "Mut","WT")
luad_egfr_del19_state <- ifelse(luad_features$Tumor_Sample_Barcode %in% substr(patients_luad_egfr_del19$Tumor_Sample_Barcode_full,1,16), "Mut","WT")
luad_egfr_mut_l858r <- ifelse(luad_features$Tumor_Sample_Barcode %in% substr(patients_luad_l858r_mut$Tumor_Sample_Barcode_full,1,16), "Yes","No")
luad_tp53_state <- ifelse(luad_features$Tumor_Sample_Barcode %in% substr(patients_luad_tp53_mut$Tumor_Sample_Barcode_full,1,16), "Mut","WT")
luad_stk11_state <- ifelse(luad_features$Tumor_Sample_Barcode %in% substr(patients_luad_stk11_mut$Tumor_Sample_Barcode_full,1,16), "Mut","WT")

luad_features$EGFR <- as.vector(luad_egfr_state)
luad_features$EGFR_L858R <- as.vector(luad_egfr_mut_l858r)
luad_features$EGFR_19Del <- as.vector(luad_egfr_del19_state)
luad_features$TP53 <- as.vector(luad_tp53_state)
luad_features$STK11 <- as.vector(luad_stk11_state)

##PDL1 groups

luad_pdl1_quantiles <- quantile(luad_features$CD274, probs = c(0.33, 0.66))

luad_pdl1_groups <- ifelse(luad_features$CD274 < luad_pdl1_quantiles[1], "Low",
                           ifelse(luad_features$CD274 >= luad_pdl1_quantiles[1] & luad_features$CD274 < luad_pdl1_quantiles[2] , "Intermediate",
                                  ifelse(luad_features$CD274 >= luad_pdl1_quantiles[2], "High", NA)))

luad_features$PDL1_expression_level <- as.vector(luad_pdl1_groups)
luad_features$PDL1_expression_level <- factor(luad_features$PDL1_expression_level,levels = c("Low","Intermediate","High"))


#Heatmap TCGA-LUAD

library(ComplexHeatmap)
library(circlize)
library(viridis)

mixture_luad_fpkm_scaled <- t(scale(t(results_mixture_luad_fpkm)))
mixture_luad_fpkm_scaled[5,] <- 0 #To avoid NAs errors due to scaling

rownames(luad_features) <- luad_features$Tumor_Sample_Barcode

col_groups <- colorRampPalette(c("green","red"))(3)

anotation_luad_mixture <-HeatmapAnnotation(df=luad_features[colnames(mixture_luad_fpkm_scaled),c(11,8,12,7,13)],
                                           annotation_name_gp = gpar(fontface="bold"),gap = unit(1,"mm"),border = TRUE,which="column",
                                           col=list(TP53= c("Mut" = "red", "WT" = "gray"),
                                                    EGFR= c("Mut" = "red", "WT" = "gray"),
                                                    STK11= c("Mut" = "red", "WT" = "gray"),
                                                    TMB_Group = c("Low"=col_groups[1],"Intermediate"=col_groups[2], "High"=col_groups[3]), 
                                                    PDL1_expression_level = c("Low"=col_groups[1],"Intermediate"=col_groups[2], "High"=col_groups[3])))

anotation_luad_mixture@anno_list[["PDL1_expression_level"]]@name <- "PDL1 mRNA"
anotation_luad_mixture@anno_list[["PDL1_expression_level"]]@color_mapping@name <- "PDL1 mRNA"

anotation_luad_mixture@anno_list[["TMB_Group"]]@name <- "TMB"
anotation_luad_mixture@anno_list[["TMB_Group"]]@color_mapping@name <- "TMB"

colour_fig_final <- colorRamp2(seq(-2,2,by=0.02),viridis(length(seq(-2,2,by=0.02))))

heatmap_luad <-Heatmap(mixture_luad_fpkm_scaled,name="Z-score",col=colour_fig_final,
                   border="black",show_column_names = FALSE,show_column_dend = FALSE,show_row_names=TRUE, show_row_dend = TRUE,
                   row_names_gp = gpar(fontsize=12), show_parent_dend_line = FALSE,row_gap = unit(0,"mm"),show_heatmap_legend = TRUE,
                   cluster_rows = TRUE,cluster_columns = TRUE,column_gap = unit(1,"mm"),column_title_gp = gpar(fontsize=12),
                   top_annotation = anotation_luad_mixture,
                   clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2")#,heatmap_legend_param = list(direction = "vertical"))
draw(heatmap_luad, merge_legend = TRUE) ##Export PDF width=10, height=6



#Barplots
library(ggpubr)

data_mixture_luad_fpkm <- t(results_mixture_luad_fpkm)
data_mixture_luad_fpkm <- merge(data_mixture_luad_fpkm, luad_features,by.x="row.names",by.y="row.names")
colnames(data_mixture_luad_fpkm)[2:23] <- gsub(" ",".", colnames(data_mixture_luad_fpkm)[2:23])

##TMB

###CD8

my_comparisons <- list( c("Low", "Intermediate"), c("Low", "High"))

p <- ggbarplot(filter(data_mixture_luad_fpkm, TMB_Group != "NA"), x = "TMB_Group", y = "T.cells.CD8", title = "TMB",
               color = "TMB_Group", line.color = "gray", line.size = 1, xlab = F,ylab = "CD8 T cells proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

###M1

my_comparisons <- list( c("Low", "Intermediate"), c("Low", "High") )

p <- ggbarplot(filter(data_mixture_luad_fpkm, TMB_Group != "NA"), x = "TMB_Group", y = "Macrophages.M1", title = "TMB",
               color = "TMB_Group", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M1 proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")


###M0

my_comparisons <- list(c("Low", "High"), c("Intermediate", "High"))

p <- ggbarplot(filter(data_mixture_luad_fpkm, TMB_Group != "NA"), x = "TMB_Group", y = "Macrophages.M0", title = "TMB",
               color = "TMB_Group", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M0 proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

###Monocytes

my_comparisons <- list(c("Low", "High"), c("Intermediate", "High") )

p <- ggbarplot(filter(data_mixture_luad_fpkm, TMB_Group != "NA"), x = "TMB_Group", y = "Monocytes", title = "TMB",
               color = "TMB_Group", line.color = "gray", line.size = 1, xlab = F,ylab = "Monocytes proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

###Dendritic activated

my_comparisons <- list( c("Low", "High"), c("Intermediate", "High") )
p <- ggbarplot(filter(data_mixture_luad_fpkm, TMB_Group != "NA"), x = "TMB_Group", y = "Dendritic.cells.activated", title = "TMB",
               color = "TMB_Group", line.color = "gray", line.size = 1, xlab = F,ylab = "Dendritic cells activated proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

###Dendritic resting

p <- ggbarplot(filter(data_mixture_luad_fpkm, TMB_Group != "NA"), x = "TMB_Group", y = "Dendritic.cells.resting", title = "TMB",
               color = "TMB_Group", line.color = "gray", line.size = 1, xlab = F,ylab = "Dendritic cells resting proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")



##EGFR

###CD8

my_comparisons <- list( c("WT", "Mut"))

p <- ggbarplot(data_mixture_luad_fpkm, x = "EGFR", y = "T.cells.CD8", title = "EGFR",
               color = "EGFR", line.color = "gray", line.size = 1, xlab = F,ylab = "CD8 T cells proportion",
               palette = c("gray", "red"),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

###Monocytes

p <- ggbarplot(data_mixture_luad_fpkm, x = "EGFR", y = "Monocytes", title = "EGFR",
               color = "EGFR", line.color = "gray", line.size = 1, xlab = F,ylab = "Monocytes proportion",
               palette = c("gray", "red"),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

###Dendritic resting

my_comparisons <- list( c("WT", "Mut"))

p <- ggbarplot(data_mixture_luad_fpkm, x = "EGFR", y = "Dendritic.cells.resting", title = "EGFR",
               color = "EGFR", line.color = "gray", line.size = 1, xlab = F,ylab = "Dendritic cells resting proportion",
               palette = c("gray", "red"),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")


##TP53

###CD8

my_comparisons <- list( c("WT", "Mut"))

p <- ggbarplot(data_mixture_luad_fpkm, x = "TP53", y = "T.cells.CD8", title = "TP53",
               color = "TP53", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD8 proportion",
               palette = c("gray", "red"),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

###M1

p <- ggbarplot(data_mixture_luad_fpkm, x = "TP53", y = "Macrophages.M1", title = "TP53",
               color = "TP53", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M1 proportion",
               palette = c("gray", "red"),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

###M0

p <- ggbarplot(data_mixture_luad_fpkm, x = "TP53", y = "Macrophages.M0", title = "TP53",
               color = "TP53", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M0 proportion",
               palette = c("gray", "red"),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

###Dendritic resting

p <- ggbarplot(data_mixture_luad_fpkm, x = "TP53", y = "Dendritic.cells.resting", title = "TP53",
               color = "TP53", line.color = "gray", line.size = 1, xlab = F,ylab = "Dendritic cells resting proportion",
               palette = c("gray", "red"),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

##PDL1 Expression

###CD8

my_comparisons <- list( c("Low", "Intermediate"), c("Low", "High"), c("Intermediate", "High") )
p <- ggbarplot(filter(data_mixture_luad_fpkm, PDL1_expression_level != "NA"), x = "PDL1_expression_level", y = "T.cells.CD8", title = "PDL1 mRNA",
               color = "PDL1_expression_level", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD8 proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

###M1

p <- ggbarplot(filter(data_mixture_luad_fpkm, PDL1_expression_level != "NA"), x = "PDL1_expression_level", y = "Macrophages.M1", title = "PDL1 mRNA",
               color = "PDL1_expression_level", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M1 proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

###Dendritic activated

my_comparisons <- list( c("Low", "High"))
p <- ggbarplot(filter(data_mixture_luad_fpkm, PDL1_expression_level != "NA"), x = "PDL1_expression_level", y = "Dendritic.cells.activated", title = "PDL1 mRNA",
               color = "PDL1_expression_level", line.color = "gray", line.size = 1, xlab = F,ylab = "Dendritic cells activated proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

##Combination of mutations

###Prepare groups

luad_egfr_stk1_tp53_comb <- ifelse(data_mixture_luad_fpkm$EGFR == "Mut" & data_mixture_luad_fpkm$TP53 != "Mut" & data_mixture_luad_fpkm$STK11 != "Mut" , "EGFR",
                                   ifelse(data_mixture_luad_fpkm$EGFR != "Mut" & data_mixture_luad_fpkm$TP53 == "Mut" & data_mixture_luad_fpkm$STK11 != "Mut", "TP53",
                                          ifelse(data_mixture_luad_fpkm$EGFR != "Mut" & data_mixture_luad_fpkm$TP53 != "Mut" & data_mixture_luad_fpkm$STK11 == "Mut", "STK11",
                                                 ifelse(data_mixture_luad_fpkm$EGFR != "Mut" & data_mixture_luad_fpkm$TP53 == "Mut" & data_mixture_luad_fpkm$STK11 == "Mut", "TP53-STK11",
                                                        ifelse(data_mixture_luad_fpkm$EGFR == "Mut" & data_mixture_luad_fpkm$TP53 != "Mut" & data_mixture_luad_fpkm$STK11 == "Mut", "EGFR-STK11",
                                                               ifelse(data_mixture_luad_fpkm$EGFR == "Mut" & data_mixture_luad_fpkm$TP53 == "Mut" & data_mixture_luad_fpkm$STK11 != "Mut", "TP53-EGFR",
                                                                      ifelse(data_mixture_luad_fpkm$EGFR == "Mut" & data_mixture_luad_fpkm$TP53 == "Mut" & data_mixture_luad_fpkm$STK11 == "Mut", "EGFR-TP53-STK11", 
                                                                             ifelse(data_mixture_luad_fpkm$EGFR != "Mut" & data_mixture_luad_fpkm$TP53 != "Mut" & data_mixture_luad_fpkm$STK11 != "Mut", "WT",NA))))))))

luad_egfr_stk1_tp53_comb <- factor(luad_egfr_stk1_tp53_comb, levels = c("WT","EGFR",
                                                                        "STK11","TP53","TP53-EGFR","TP53-STK11","EGFR-TP53-STK11"))

luad_egfr_del_l858r <- ifelse(data_mixture_luad_fpkm$EGFR_19Del == "Mut" & data_mixture_luad_fpkm$EGFR_L858R != "Yes" , "Del19",
                              ifelse(data_mixture_luad_fpkm$EGFR_19Del != "Mut" & data_mixture_luad_fpkm$EGFR_L858R == "Yes", "L858R",
                                     ifelse(data_mixture_luad_fpkm$EGFR_19Del != "Mut" & data_mixture_luad_fpkm$EGFR_L858R != "Yes" & data_mixture_luad_fpkm$EGFR == "Mut", "Other",
                                            ifelse(data_mixture_luad_fpkm$EGFR != "Mut", "WT", NA))))

luad_egfr_del_l858r <- factor(luad_egfr_del_l858r, levels = c("WT","Other","L858R","Del19"))

data_mixture_luad_fpkm$EGFR_STK11_TP53 <- luad_egfr_stk1_tp53_comb
data_mixture_luad_fpkm$EGFR_DEL19_L858R <- luad_egfr_del_l858r

###Combinations with EGFR, STK11, TP53

###CD8

my_comparisions_egfr_stk11_tp53_comb <- list(c("TP53", "WT"),c("WT","TP53-EGFR"), c("EGFR", "TP53"), c("EGFR", "TP53-STK11"),
                                             c("TP53", "TP53-EGFR"),
                                             c("TP53-STK11", "TP53-EGFR"))

p <- ggbarplot(filter(data_mixture_luad_fpkm, EGFR_STK11_TP53!="EGFR-TP53-STK11"), x = "EGFR_STK11_TP53", y = "T.cells.CD8", title = "Mutations",
               color = "EGFR_STK11_TP53", line.color = "gray", line.size = 1, xlab = F,ylab = "CD8 T cells proportion",
               short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisions_egfr_stk11_tp53_comb, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

###M1

my_comparisions_egfr_stk11_tp53_comb <- list( c("TP53", "WT"),c("EGFR", "TP53"),c("EGFR","TP53-EGFR"),
                                              c("EGFR", "TP53-STK11"), c("TP53","STK11"))

p <- ggbarplot(filter(data_mixture_luad_fpkm, EGFR_STK11_TP53!="EGFR-TP53-STK11"), x = "EGFR_STK11_TP53", y = "Macrophages.M1", title = "Mutations",
               color = "EGFR_STK11_TP53", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M1 proportion",
               short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisions_egfr_stk11_tp53_comb, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

###Monocytes

my_comparisions_egfr_stk11_tp53_comb <- list(c("EGFR","WT"),c("EGFR","STK11"),c("EGFR", "TP53"), c("EGFR", "TP53-STK11"))

p <- ggbarplot(filter(data_mixture_luad_fpkm, EGFR_STK11_TP53!="EGFR-TP53-STK11"), x = "EGFR_STK11_TP53", y = "Monocytes", title = "Mutations",
               color = "EGFR_STK11_TP53", line.color = "gray", line.size = 1, xlab = F,ylab = "Monocytes proportion",
               short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisions_egfr_stk11_tp53_comb, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")


###EGFR deletions in exon 19 and L858R mutations

###CD8

my_comparisions_egfr_del_l858r <- list(c("Del19", "WT"),c("L858R","WT"))

p <- ggbarplot(data_mixture_luad_fpkm, x = "EGFR_DEL19_L858R", y = "T.cells.CD8", title = "EGFR mutations",
               color = "EGFR_DEL19_L858R", line.color = "gray", line.size = 1, xlab = F,ylab = "CD8 T cells proportion",
               short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisions_egfr_del_l858r, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8,hide.ns = F)
ggpar(p, legend="none")

###M1

my_comparisions_egfr_del_l858r <- list(c("Del19", "WT"))

p <- ggbarplot(data_mixture_luad_fpkm, x = "EGFR_DEL19_L858R", y = "Macrophages.M1", title = "EGFR mutations",
               color = "EGFR_DEL19_L858R", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M1 proportion",
               short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisions_egfr_del_l858r, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8,hide.ns = T)
ggpar(p, legend="none")


###Monocytes

my_comparisions_egfr_del_l858r <- list( c("Del19", "WT"),c("L858R","WT"))

p <- ggbarplot(data_mixture_luad_fpkm, x = "EGFR_DEL19_L858R", y = "Monocytes", title = "EGFR mutations",
               color = "EGFR_DEL19_L858R", line.color = "gray", line.size = 1, xlab = F,ylab = "Monocytes proportion",
               short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisions_egfr_del_l858r, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8,hide.ns = T)
ggpar(p, legend="none")


#Tregs

colnames(data_mixture_luad_fpkm)[10] <- "T.cells.regulatory"

my_comparisions_egfr_del_l858r <- list( c("Del19", "WT"))

p <- ggbarplot(data_mixture_luad_fpkm, x = "EGFR_DEL19_L858R", y = "T.cells.regulatory", title = "EGFR mutations",
               color = "EGFR_DEL19_L858R", line.color = "gray", line.size = 1, xlab = F,ylab = "T Regs proportion",
               short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisions_egfr_del_l858r, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8,hide.ns = T)
ggpar(p, legend="none")


###T cells CD4 memory activated

my_comparisions_egfr_del_l858r <- list( c("Del19", "Other"), c("Del19", "WT"),c("L858R","WT"))

p <- ggbarplot(data_mixture_luad_fpkm, x = "EGFR_DEL19_L858R", y = "T.cells.CD4.memory.activated", title = "EGFR mutations",
               color = "EGFR_DEL19_L858R", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD4 memory activated proportion",
               short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisions_egfr_del_l858r, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8,hide.ns = T)
ggpar(p, legend="none")


###T cells CD4 memory resting

my_comparisions_egfr_del_l858r <- list(c("Del19", "Other"), c("Del19", "WT"),c("L858R","Other"),c("L858R","WT"))
p <- ggbarplot(data_mixture_luad_fpkm, x = "EGFR_DEL19_L858R", y = "T.cells.CD4.memory.resting", title = "EGFR mutations",
               color = "EGFR_DEL19_L858R", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD4 memory resting proportion",
               short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisions_egfr_del_l858r, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8,hide.ns = T)
ggpar(p, legend="none")




#TCGA-SKCM analysis

##Filter normal samples and replicates

tumor_samples_skcm <- colnames(assay_data_melanoma_df)[substr(colnames(assay_data_melanoma_df),14,15) !=11]

assay_data_melanoma_df <- assay_data_melanoma_df[,tumor_samples_skcm]


##MIXTURE
sgm <- assay_data_melanoma_df
source("../../Utils/RUN_MIXTURE.R") #Check RUN_MIXTURE.R location - This run MIXTURE
#mix.test object is generated after run
results_mixture_skcm_fpkm <-  mix.test$Subjects$MIXprop #relative proportions
results_mixture_skcm_fpkm <- t(results_mixture_skcm_fpkm)
results_mixture_skcm_fpkm <- as.data.frame(results_mixture_skcm_fpkm)

##ITH data

data_skcm_ith<- read.table("./ith.txt",sep = "\t",header = T) #ITH data from Wolf et al 2019

skcm_features <- data.frame(Barcode = colnames(results_mixture_skcm_fpkm),Patient=substr(colnames(results_mixture_skcm_fpkm),1,12))
skcm_features <- merge(skcm_features,data_skcm_ith,by.x="Patient",by.y="Patient")

skcm_features$clones <- as.factor(skcm_features$clones)
names(skcm_features)[8] <- "Cytolytic Score"
clones_3_groups_skcm <- skcm_features$clones
clones_3_groups_skcm <- ifelse(clones_3_groups_skcm == 1 | clones_3_groups_skcm == 2, "Low",
                               ifelse(clones_3_groups_skcm == 3 | clones_3_groups_skcm == 4, "Intermediate",
                                      ifelse(clones_3_groups_skcm == 5 | clones_3_groups_skcm == 6, "High", NA)))
skcm_features$ITH <- as.vector(clones_3_groups_skcm)
skcm_features$ITH <- factor(skcm_features$ITH,levels = c("Low","Intermediate","High"))
rownames(skcm_features) <- skcm_features$Barcode

##Heatmap
data_mixture_skcm_scaled <- t(scale(t(results_mixture_skcm_fpkm[,rownames(skcm_features)])))

anotation_skcm_mixture <-HeatmapAnnotation(df=skcm_features[,c(13,8)],
                                           annotation_name_gp = gpar(fontface="bold"),gap = unit(1,"mm"),border = TRUE,which="column",
                                           col=list(ITH = c("Low"=col_groups[1],"Intermediate"=col_groups[2], "High"=col_groups[3]),
                                           "Cytolytic Score" = colorRamp2(c(0,2000,4000,6000), c("#FFFFFF", "#F4C2E3", "#E483C8", "#DB63BB"))))

heatmap_skcm<-Heatmap(data_mixture_skcm_scaled,name="Z-score",col=colour_fig_final,
                   border="black",show_column_names = FALSE,show_column_dend = FALSE,show_row_names=TRUE, show_row_dend = TRUE,
                   row_names_gp = gpar(fontsize=12), show_parent_dend_line = FALSE,row_gap = unit(0,"mm"),show_heatmap_legend = TRUE,
                   cluster_rows = TRUE,cluster_columns = TRUE,column_gap = unit(1,"mm"),column_title_gp = gpar(fontsize=12),
                   top_annotation = anotation_skcm_mixture,
                   clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2")
draw(heatmap_skcm, merge_legend = TRUE) ##Export PDF width=10, height=6



data_skcm_plots <-as.data.frame(t(results_mixture_skcm_fpkm[,rownames(skcm_features)]))
data_skcm_plots$ITH <- skcm_features$ITH
data_skcm_plots$`Cytolytic Score` <- skcm_features$`Cytolytic Score`
colnames(data_skcm_plots)[1:22] <- gsub(" ",".", colnames(data_skcm_plots)[1:22])

skcm_cyt_quantiles <- quantile(data_skcm_plots$`Cytolytic Score`, probs = c(0.33, 0.66))

skcm_cyt_groups <- ifelse(data_skcm_plots$`Cytolytic Score` < skcm_cyt_quantiles[1], "Low",
                           ifelse(data_skcm_plots$`Cytolytic Score` >= skcm_cyt_quantiles[1] & data_skcm_plots$`Cytolytic Score` < skcm_cyt_quantiles[2] , "Intermediate",
                                  ifelse(data_skcm_plots$`Cytolytic Score` >= skcm_cyt_quantiles[2], "High", NA)))

data_skcm_plots$CYT <- as.vector(skcm_cyt_groups)
data_skcm_plots$CYT <- factor(data_skcm_plots$CYT,levels = c("Low","Intermediate","High"))

##Plots

### ITH

#### CD8

my_comparisons <- list( c("Low", "Intermediate"), c("Low", "High"))

p <- ggbarplot(data_skcm_plots, x = "ITH", y = "T.cells.CD8", title = "ITH",
               color = "ITH", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD8 proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8,hide.ns = F)
ggpar(p, legend="none")

#### M2

my_comparisons <- list(c("Low", "High"))
p <- ggbarplot(data_skcm_plots, x = "ITH", y = "Macrophages.M2", title = "ITH",
               color = "ITH", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M2 proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8,hide.ns = F)
ggpar(p, legend="none")

#### T cells CD4 memory activated

my_comparisons <- list(c("Low", "High"), c("Intermediate", "High") )
p <- ggbarplot(data_skcm_plots, x = "ITH", y = "T.cells.CD4.memory.activated", title = "ITH",
               color = "ITH", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD4 memory activated proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

#### T cells CD4 memory resting

my_comparisons <- list(c("Low", "High"))
p <- ggbarplot(data_skcm_plots, x = "ITH", y = "T.cells.CD4.memory.resting", title = "ITH",
               color = "ITH", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD4 memory resting proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")



### Cytolytic Score

#### CD8

my_comparisons <- list(c("Low","Intermediate"),c("Low", "High"),c("Intermediate","High"))
p <- ggbarplot(data_skcm_plots, x = "CYT", y = "T.cells.CD8", title = "Cytolytic Score",
               color = "CYT", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD8 proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

#### M0

my_comparisons <- list(c("Low","Intermediate"),c("Low", "High"),c("Intermediate","High"))
p <- ggbarplot(data_skcm_plots, x = "CYT", y = "Macrophages.M0", title = "Cytolytic Score",
               color = "CYT", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M0 proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

#### M1

my_comparisons <- list(c("Low","Intermediate"),c("Low", "High"),c("Intermediate","High"))
p <- ggbarplot(data_skcm_plots, x = "CYT", y = "Macrophages.M1", title = "Cytolytic Score",
               color = "CYT", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M1 proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

#### M2

my_comparisons <- list(c("Low","Intermediate"),c("Low", "High"),c("Intermediate","High"))
p <- ggbarplot(data_skcm_plots, x = "CYT", y = "Macrophages.M2", title = "Cytolytic Score",
               color = "CYT", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M2 proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

#### B cells naive

my_comparisons <- list(c("Low", "High"))
p <- ggbarplot(data_skcm_plots, x = "CYT", y = "B.cells.naive", title = "Cytolytic Score",
               color = "CYT", line.color = "gray", line.size = 1, xlab = F,ylab = "B cells naive proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

#### Plasma cells

my_comparisons <- list(c("Low","Intermediate"),c("Low", "High"),c("Intermediate","High"))
p <- ggbarplot(data_skcm_plots, x = "CYT", y = "Plasma.cells", title = "Cytolytic Score",
               color = "CYT", line.color = "gray", line.size = 1, xlab = F,ylab = "Plasma cells proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

#### T cells CD4 memory activated

my_comparisons <- list(c("Low","Intermediate"),c("Low", "High"),c("Intermediate","High"))
p <- ggbarplot(data_skcm_plots, x = "CYT", y = "T.cells.CD4.memory.activated", title = "Cytolytic Score",
               color = "CYT", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD4 memory activated proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

#### T cells CD4 memory resting

my_comparisons <- list(c("Low","Intermediate"),c("Low", "High"),c("Intermediate","High"))
p <- ggbarplot(data_skcm_plots, x = "CYT", y = "T.cells.CD4.memory.resting", title = "Cytolytic Score",
               color = "CYT", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD4 memory resting proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

#### T regs

colnames(data_skcm_plots)[9] <- "T.cells.regulatory"

my_comparisons <- list(c("Low","Intermediate"),c("Low", "High"))
p <- ggbarplot(data_skcm_plots, x = "CYT", y = "T.cells.regulatory", title = "Cytolytic Score",
               color = "CYT", line.color = "gray", line.size = 1, xlab = F,ylab = "T regs proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

#### NK activated

my_comparisons <- list(c("Low","Intermediate"),c("Low", "High"))
p <- ggbarplot(data_skcm_plots, x = "CYT", y = "NK.cells.activated", title = "Cytolytic Score",
               color = "CYT", line.color = "gray", line.size = 1, xlab = F,ylab = "NK cells activated proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")


#TCGA-HNSC analysis

##Filter normal samples

tumor_samples_hnsc <- colnames(assay_data_HNSC_df)[substr(colnames(assay_data_HNSC_df),14,15) !=11]

assay_data_HNSC_df <- assay_data_HNSC_df[,tumor_samples_hnsc]

##MIXTURE
sgm <- assay_data_HNSC_df
source("../../Utils/RUN_MIXTURE.R") #Check RUN_MIXTURE.R location - This run MIXTURE
#mix.test object is generated after run
results_mixture_hnsc_fpkm <-  mix.test$Subjects$MIXprop #relative proportions
results_mixture_hnsc_fpkm <- t(results_mixture_hnsc_fpkm)
results_mixture_hnsc_fpkm <- as.data.frame(results_mixture_hnsc_fpkm)

##ITH

tcga_load("HNSC")
math_score_hnsc <- math.score(maf = tcga_hnsc_mc3)
hnsc_features <- data.frame(Sample = colnames(results_mixture_hnsc_fpkm), Tumor_Sample_Barcode = substr(colnames(results_mixture_hnsc_fpkm),1,12))
hnsc_features <- merge(hnsc_features,math_score_hnsc,by= "Tumor_Sample_Barcode", all.x=T)
rownames(hnsc_features) <- hnsc_features$Sample
hnsc_features <- hnsc_features[,-2]

hnsc_math_quantiles <- quantile(hnsc_features$MATH, probs = c(0.33, 0.66),na.rm=T)
hnsc_math_groups <- ifelse(hnsc_features$MATH < hnsc_math_quantiles[1], "Low",
                           ifelse(hnsc_features$MATH >= hnsc_math_quantiles[1] & hnsc_features$MATH < hnsc_math_quantiles[2] , "Intermediate",
                                  ifelse(hnsc_features$MATH >= hnsc_math_quantiles[2], "High", NA)))

hnsc_features$ITH <- as.vector(hnsc_math_groups)
hnsc_features$ITH <- factor(hnsc_features$ITH,levels = c("Low","Intermediate","High"))

##Plots

data_mixture_hnsc <- merge(t(results_mixture_hnsc_fpkm),hnsc_features,by="row.names")
colnames(data_mixture_hnsc)[2:23] <- gsub(" ",".",colnames(data_mixture_hnsc)[2:23])


### CD8

my_comparisons <- list( c("Low", "Intermediate"), c("Low", "High"), c("Intermediate", "High") )

p <- ggbarplot(filter(data_mixture_hnsc,ITH != "NA"), x = "ITH", y = "T.cells.CD8", title = "ITH",
               color = "ITH", line.color = "gray", line.size = 1, xlab = F,ylab = "CD8 T cells proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F,
               add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")


### M0

my_comparisons <- list( c("Low", "Intermediate"), c("Low", "High") )
p <- ggbarplot(filter(data_mixture_hnsc,ITH != "NA"), x = "ITH", y = "Macrophages.M0", title = "ITH",
               color = "ITH", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M0 proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

### M1

my_comparisons <- list(  c("Low", "High") )
p <- ggbarplot(filter(data_mixture_hnsc,ITH != "NA"), x = "ITH", y = "Macrophages.M1", title = "ITH",
               color = "ITH", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M1 proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

### CD4 memory activated

my_comparisons <- list( c("Low", "Intermediate"), c("Low", "High") )
p <- ggbarplot(filter(data_mixture_hnsc,ITH != "NA"), x = "ITH", y = "T.cells.CD4.memory.activated", title = "ITH",
               color = "ITH", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD4 memory activated proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")

### CD4 memory resting

my_comparisons <- list( c("Low", "Intermediate"), c("Low", "High"))
p <- ggbarplot(filter(data_mixture_hnsc,ITH != "NA"), x = "ITH", y = "T.cells.CD4.memory.resting", title = "ITH",
               color = "ITH", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD4 memory resting proportion",
               palette = c(col_groups[1],col_groups[2],col_groups[3]),short.panel.labs = T,combine=F, add = c("mean_sd","jitter"))+theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", paired = FALSE,size=8)
ggpar(p, legend="none")



#TCGA-COAD analysis

##Filter normal samples

tumor_samples_coad <- colnames(assay_data_COAD_df)[substr(colnames(assay_data_COAD_df),14,15) !=11]
tumor_samples_coad_replicates <- colnames(assay_data_COAD_df)[duplicated(substr(colnames(assay_data_COAD_df),1,16))]

assay_data_COAD_df <- assay_data_COAD_df[,tumor_samples_coad]
assay_data_COAD_df <- assay_data_COAD_df[,!colnames(assay_data_COAD_df) %in% tumor_samples_coad_replicates]
assay_data_COAD_counts_df <- assay_data_COAD_counts_df[,tumor_samples_coad]
assay_data_COAD_counts_df <- assay_data_COAD_counts_df[,!colnames(assay_data_COAD_counts_df) %in% tumor_samples_coad_replicates]

##MIXTURE
sgm <- assay_data_COAD_df
source("../../Utils/RUN_MIXTURE.R") #Check RUN_MIXTURE.R location - This run MIXTURE
#mix.test object is generated after run
results_mixture_coad_fpkm <-  mix.test$Subjects$MIXprop #relative proportions
results_mixture_coad_fpkm <- t(results_mixture_coad_fpkm)
results_mixture_coad_fpkm <- as.data.frame(results_mixture_coad_fpkm)

##MSI STATUS

###MSI data downloaded from: http://krishna.gs.washington.edu/content/members/hauser/mosaic/processed_exome_microsatellite_instability_calls/

load("./batch1_results_post_review_msi_predicted_063016.robj")
load("./batch2_results_post_review_msi_predicted_070516.robj")
sampledat2n_1 <- sampledat2n
load("./batch2_results_post_review_msi_predicted_070616.robj")
load("./batch3_results_post_review_msi_predicted_070716.robj")
load("./batch4_results_post_review_msi_predicted_070516.robj")
load("./batch5_results_post_review_msi_predicted_070516.robj")
load("./batch6_results_post_review_msi_predicted_070516.robj")

msi_status <- rbind(sampledat1n,sampledat2n_1,sampledat2n,sampledat3n,sampledat4n,sampledat5n,sampledat6n)

## Sample data
features_coad <- as.data.frame(t(as.data.frame(colData_COAD)))
colnames(features_coad) <- unlist(features_coad["sample",])
features_coad <- features_coad[,tumor_samples_coad]
features_coad <- features_coad[,!colnames(features_coad) %in% tumor_samples_coad_replicates]
features_coad <- as.data.frame(t(features_coad))
features_coad <- features_coad[,c("sample","patient","paper_MSI_status")]
features_coad <- data.frame(sample = unlist(features_coad$sample), patient= unlist(features_coad$patient), paper_MSI_status=unlist(features_coad$paper_MSI_status))
rownames(features_coad) <- NULL
features_coad$paper_MSI_status <- as.character(features_coad$paper_MSI_status)
features_coad$sample <- as.character(features_coad$sample)
features_coad$patient <- as.character(features_coad$patient)

###Merge MSI data
msi_status_coad <- filter(msi_status, sample_name %in% features_coad$patient)
msi_status_coad$sample_name <- as.character(msi_status_coad$sample_name)
features_coad <- merge(features_coad,msi_status_coad,by.x="patient",by.y="sample_name",all.x=T)
features_coad <- mutate(features_coad, MSI = ifelse(!is.na(features_coad$msi),features_coad$msi,
                                                    ifelse(!is.na(features_coad$paper_MSI_status),features_coad$paper_MSI_status,NA)))
features_coad$MSI[features_coad$MSI == "MSI-L"] = "MSS"
patient_coad_msi_h <- filter(features_coad, MSI == "MSI-H") %>% dplyr::select(sample) %>% unlist()
patient_coad_mss <- filter(features_coad, MSI == "MSS") %>% dplyr::select(sample) %>% unlist()

##CYT

###Cytolitic score

####GZMA and PRF1 expression

coad_expression_tmm <- DGEList(counts=assay_data_COAD_counts_df)
coad_expression_tmm <- calcNormFactors(coad_expression_tmm, method = "TMM")
tmm_coad <- cpm(coad_expression_tmm)
tmm_coad <- tmm_coad[c("GZMA","PRF1"),]
tmm_coad <- t(tmm_coad)
tmm_coad <- as.data.frame(tmm_coad)
tmm_coad$Tumor_Sample_Barcode <- rownames(tmm_coad)
rownames(tmm_coad) <- NULL
tmm_coad <- tmm_coad[,c(3,1,2)]
#tmm_coad <- mutate(tmm_coad, short_sample_code = substr(tmm_coad$Tumor_Sample_Barcode,1,15))
tmm_coad <- mutate(tmm_coad, `Cytolytic Score`= sqrt(GZMA*PRF1))

features_coad$MSI <- factor(features_coad$MSI,levels = c("MSS","MSI-H"))

features_coad <- merge(features_coad,tmm_coad,by.x="sample",by.y="Tumor_Sample_Barcode", all.x=T)

##Heatmap

col_cyt_coad = colorRamp2(c(0, 20, 40,60), c("#FFFFFFFF", "#F4C2E3FF", "#E483C8FF", "#DB63BBFF"))

features_coad_noNA <- filter(features_coad,MSI != "NA")
rownames(features_coad_noNA) <- features_coad_noNA$sample


results_mixture_COAD_scaled <-  t(scale(t(results_mixture_coad_fpkm)))
results_mixture_COAD_scaled[5,] <- 0
colnames(results_mixture_COAD_scaled) <- substr(colnames(results_mixture_COAD_scaled),1,16)
results_mixture_COAD_scaled <- results_mixture_COAD_scaled[,rownames(features_coad_noNA)]


ha_coad = HeatmapAnnotation(MSI = features_coad_noNA$MSI,`Cytolytic Score` = features_coad_noNA$`Cytolytic Score`,
                       col=list(MSI=c("MSS" = "blue", "MSI-H" = "red"),`Cytolytic Score`=col_cyt_coad),
                       annotation_name_gp = gpar(fontface="bold"),gap = unit(1,"mm"),border = TRUE,which="column")

plot_coad_heatmap<-Heatmap(results_mixture_COAD_scaled,name="Z-score",col=colour_fig_final,
                           border="black",show_column_names = FALSE,show_column_dend = F,show_row_names=TRUE, show_row_dend = TRUE,
                           row_names_gp = gpar(fontsize=12), show_parent_dend_line = FALSE,row_gap = unit(0,"mm"),show_heatmap_legend = TRUE,
                           cluster_rows = TRUE,cluster_columns = TRUE,column_gap = unit(1,"mm"),column_title_gp = gpar(fontsize=12),
                           top_annotation = ha_coad, column_split = features_coad_noNA$MSI,
                           clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2")
draw(plot_coad_heatmap, merge_legend = TRUE) ##Export PDF width=10, height=6

##Plots
colnames(results_mixture_coad_fpkm) <- substr(colnames(results_mixture_coad_fpkm),1,16)
data_mixture_coad <- merge(t(results_mixture_coad_fpkm),features_coad_noNA,by="row.names")
colnames(data_mixture_coad)[2:23] <- gsub(" ","_",colnames(data_mixture_coad)[2:23])

### CD8
my_comparisons <- list( c("MSS", "MSI-H"))
p <- ggbarplot(data_mixture_coad, x = "MSI", y = "T_cells_CD8", title = "MSI",
               color = "MSI", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD8 proportion",
               palette = c("blue","red"),short.panel.labs = T,combine=F, add = c("mean_se","jitter")) + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox", paired = FALSE,size=8)
ggpar(p, legend="none")

### Macrophages M1

p <- ggbarplot(data_mixture_coad, x = "MSI", y = "Macrophages_M1", title = "MSI",
               color = "MSI", line.color = "gray", line.size = 1, xlab = F,ylab = "Macrophages M1 proportion",
               palette = c("blue","red"),short.panel.labs = T,combine=F, add = c("mean_se","jitter")) + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox", paired = FALSE,size=8)
ggpar(p, legend="none")

### CD4 memory activated

p <- ggbarplot(data_mixture_coad, x = "MSI", y = "T_cells_CD4_memory_activated", title = "MSI",
               color = "MSI", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD4 memory activated proportion",
               palette = c("blue","red"),short.panel.labs = T,combine=F, add = c("mean_se","jitter")) + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox", paired = FALSE,size=8)
ggpar(p, legend="none")


### CD4 memory resting

p <- ggbarplot(data_mixture_coad, x = "MSI", y = "T_cells_CD4_memory_resting", title = "MSI",
               color = "MSI", line.color = "gray", line.size = 1, xlab = F,ylab = "T cells CD4 memory resting proportion",
               palette = c("blue","red"),short.panel.labs = T,combine=F, add = c("mean_se","jitter")) + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox", paired = FALSE,size=8)
ggpar(p, legend="none")


### NK activated

p <- ggbarplot(data_mixture_coad, x = "MSI", y = "NK_cells_activated", title = "MSI",
               color = "MSI", line.color = "gray", line.size = 1, xlab = F,ylab = "NK cells activated proportion",
               palette = c("blue","red"),short.panel.labs = T,combine=F, add = c("mean_se","jitter")) + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox", paired = FALSE,size=8)
ggpar(p, legend="none")


### NK resting

p <- ggbarplot(data_mixture_coad, x = "MSI", y = "NK_cells_resting", title = "MSI",
               color = "MSI", line.color = "gray", line.size = 1, xlab = F,ylab = "NK cells resting proportion",
               palette = c("blue","red"),short.panel.labs = T,combine=F, add = c("mean_se","jitter")) + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox", paired = FALSE,size=8)
ggpar(p, legend="none")


### Plasma cells

p <- ggbarplot(data_mixture_coad, x = "MSI", y = "Plasma_cells", title = "MSI",
               color = "MSI", line.color = "gray", line.size = 1, xlab = F,ylab = "Plasma cells proportion",
               palette = c("blue","red"),short.panel.labs = T,combine=F, add = c("mean_se","jitter")) + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox", paired = FALSE,size=8)
ggpar(p, legend="none")


### Neutrophils

p <- ggbarplot(data_mixture_coad, x = "MSI", y = "Neutrophils", title = "MSI",
               color = "MSI", line.color = "gray", line.size = 1, xlab = F,ylab = "Neutrophils proportion",
               palette = c("blue","red"),short.panel.labs = T,combine=F, add = c("mean_se","jitter")) + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox", paired = FALSE,size=8)
ggpar(p, legend="none")


### Mast cells resting

p <- ggbarplot(data_mixture_coad, x = "MSI", y = "Mast_cells_resting", title = "MSI",
               color = "MSI", line.color = "gray", line.size = 1, xlab = F,ylab = "Mast cells resting proportion",
               palette = c("blue","red"),short.panel.labs = T,combine=F, add = c("mean_se","jitter")) + theme(text = element_text(size = 20),plot.title = element_text(hjust = 0.5))  +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox", paired = FALSE,size=8)
ggpar(p, legend="none")
