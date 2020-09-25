#Supplementary Figure_3: Performance of genset enrichment analysis methods on pure cell lines data
library(ggplot2)
library(ggpubr)

rm(list=ls())

# Generate plotting function
gplot <- function(df, title){
  ggplot(df, aes(y=Score,x=cells)) + geom_boxplot() +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, size=8,hjust = 1),
          axis.title.x = element_blank()) +
    labs(title = title, y = paste(unique(df$Method),"Score"))
}

## xCell results

# Load data
xcell.estimates <- read.table("Figure_2/Supplementary_Figure3/xCell_estimates.txt",h=T, sep="\t",row.names = "X")
xcell.pvals <- read.table("Figure_2/Supplementary_Figure3/xCell_pvalues.txt",h=T, sep="\t",row.names = "X")
xcell.raw <- read.table("Figure_2/Supplementary_Figure3/xCell_RAW.txt", h=T, sep="\t",row.names = "X")

# Extract xCell beta and p-values
betas <- t(xcell.estimates)
colnames(betas)<-rownames(xcell.estimates)

pvalues <- t(xcell.pvals)
colnames(pvalues) <- rownames(xcell.pvals)

# Set non-significant gene sets to NA
betas_ns <- betas
betas_ns[pvalues > 0.05] <- NA

# Take a look at distributions
boxplot(betas_ns)
dim(betas_ns)

# Generate subgroups of cell types
Lymphoids <- c("B-cells","Memory B-cells","naive B-cells","pro B-cells","Class-switched memory B-cells",
               "CD4+ T-cells","CD4+ Tcm","CD4+ Tem","CD4+ memory T-cells","CD4+ naive T-cells",
               "CD8+ T-cells","CD8+ Tcm","CD8+ Tem","CD8+ naive T-cells","Tgd cells","Th1 cells",
               "Th2 cells","Tregs","NK cells","NKT","Plasma cells")
Myeloids <- c("Basophils","Eosinophils","DC","cDC","iDC","pDC","aDC",
              "Monocytes","Macrophages","Macrophages M1","Macrophages M2",
              "Mast cells","Neutrophils")
StemCells <- c("MSC","HSC","CLP","CMP","GMP","MEP","MPP","Megakaryocytes","Erythrocytes","Platelets")
StromaCells <- c("Adipocytes","Fibroblasts","Pericytes","Endothelial cells",
                 "Chondrocytes","Smooth muscle","Osteoblast","Myocytes","Skeletal muscle",
                 "ly Endothelial cells","mv Endothelial cells","Preadipocytes")
Others <- c("Astrocytes","Mesangial cells","Hepatocytes","Epithelial cells",
           "Keratinocytes","Melanocytes","Sebocytes","Neurons")
Scores <- c("ImmuneScore","StromaScore","MicroenvironmentScore")

# Extract betas for lymphoid and myeloid cells
xCell_myeloids_lymphoids <- betas_ns[, c(Lymphoids, Myeloids)]
dim(xCell_myeloids_lymphoids)
xCell_myeloids_lymphoids_df <- data.frame(Score= c(xCell_myeloids_lymphoids),cells=rep(c(Lymphoids, Myeloids),each=nrow(xCell_myeloids_lymphoids)))
xCell_myeloids_lymphoids_df$Method <- "xCell"

# Generate plot
xCell.plot <- gplot(df = xCell_myeloids_lymphoids_df,title = "xCell")
xCell.plot
dev.off()

## IMSIG results

# Load data
IMSIG.estimates <- read.table("Figure_2/Supplementary_Figure3/ImSig_estimates.txt", h=T)
colnames(IMSIG.estimates)

# Take a look at distributions
boxplot(IMSIG.estimates[,c(1,3,4,5,6,7,9)])
summary(rowSums(IMSIG.estimates[,c(1,3,4,5,6,7,9)]>0,na.rm=T))
colnames(IMSIG.estimates)

# Define a subset of immune cell types
IMSIG_immunes <- c("B.cells","Macrophages","Monocytes","Neutrophils","NK.cells","Plasma.cells","T.cells")
IMSIG.estimates.immune <- IMSIG.estimates[,IMSIG_immunes]

# Generate data frame
IMSIG_immune_df <- data.frame(Score = c(data.matrix(IMSIG.estimates.immune)), cells = rep(colnames(IMSIG.estimates.immune),each=nrow(IMSIG.estimates.immune)))
IMSIG_immune_df$Method <- "ImSig"
IMSIG_immune_df$cells <- stringr::str_replace_all(IMSIG_immune_df$cells,"\\."," ")

# Generate plot
IMSIG.plot <- gplot(IMSIG_immune_df,"ImSig") + ylim(0, max(IMSIG_immune_df$Score))
IMSIG.plot
dev.off()

## MCPcounter results

# Load data
mcp.estimates <- t(read.table("Figure_2/Supplementary_Figure3/MCP_estimates.txt", h=T))
colnames(mcp.estimates)
summary(rowSums(mcp.estimates>0))

# Define a subset of immune cell types
mcp_immunes <- c("T cells","CD8 T cells","Cytotoxic lymphocytes","NK cells","B lineage","Monocytic lineage",
                "Myeloid dendritic cells","Neutrophils")

# Generate data frame
mcp_estimates_immune_df <- data.frame(Score = c(data.matrix(mcp.estimates[, mcp_immunes])), cells = rep(colnames(mcp.estimates[, mcp_immunes]),each=nrow(mcp.estimates)))
mcp_estimates_immune_df$Method <- "MCP_counter"

#Generate plot
mcp.plot <- gplot(mcp_estimates_immune_df,"MCPcounter") + ylim(0, max(mcp_estimates_immune_df$Score))
mcp.plot
dev.off()

# Save plots in PDF 16x5
ggsave(ggarrange(plotlist = list(xCell.plot, IMSIG.plot, mcp.plot),ncol = 3,widths= c(2,1,1),align = "h"),filename = "Figure_2/Supplementary_Figure3/Supplementary_Figure3.pdf",device = "pdf",width = 16,height = 5)







