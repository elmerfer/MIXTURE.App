#Figure_2e-h: Performance in pure cell lines data set, cell types estimation (e-f) and proportions (g-h)
rm(list=ls())
library(gridExtra)
library(ggplot2)

#Load matrix with MIXTURE output using LM22 signature on CCLE data
df.lm22 <- readRDS("Figure_2/Figure_2e-h/CellLinesPlotLM22.Rds")
#Load matrix with MIXTURE output using TIl10 signature on CCLE data
df.til10 <- readRDS("Figure_2/Figure_2e-h/CellLinesPlotTIL10.Rds")

#Binding and formating of both matrices
df.t <- rbind(df.lm22, df.til10)
colnames(df.t)[2] <- c("Total")
df.t$Method <- gsub("QUANTISEQ","quanTIseq",df.t$Method)
df.t$Method <- factor(df.t$Method,levels=c("ABBAS","ABIS","quanTIseq","CIBERSORT","MIXTURE"))

#Plot estimated cell types
PLOT1 <- ggplot(df.t[c(1:5090),], aes(Method, NumCels, fill = Method)) +  geom_boxplot(outlier.size = 0.0) + geom_jitter(height = 0, width = 0.2, size = 0.3) +
  theme_classic() + ylab(label="Estimated number of cell types") + theme(legend.position = "none",plot.title =element_text(hjust = 0.5,size = 30),axis.text.y.left = element_text(size = 28),
                                                                         axis.title.y.left = element_text(size = 30),axis.text.x = element_text(size = 28,angle=45,vjust = 0.4,hjust = 0.5),
                                                                         axis.ticks.length = unit(0.4,"cm"),axis.ticks = element_line(size = 1)) + xlab(label="") + ggtitle(label="LM22")
PLOT2 <- ggplot(df.t[c(5091:10180),], aes(Method, NumCels, fill = Method)) +  geom_boxplot(outlier.size = 0.0) + geom_jitter(height = 0, width = 0.2, size = 0.3) +
  theme_classic() + ylab(label="Estimated number of cell types") + theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 30),axis.text.y.left = element_text(size = 28),
                                                                         axis.title.y.left = element_text(size = 30),axis.text.x = element_text(size = 28,angle=45,vjust = 0.4,hjust = 0.5),
                                                                         axis.ticks.length = unit(0.4,"cm"),axis.ticks = element_line(size = 1)) + xlab(label="") + ggtitle(label = "TIL10")

#Plot estimated proportions
PLOT3 <- ggplot(df.t[c(1:5090),], aes(Method,Total, fill = Method)) + geom_violin(scale = "width") + geom_boxplot(width = 0.1,outlier.size = 0.5) + ggtitle(label="LM22") +
  theme_classic() + ylab(label="Total Unnormalized Immune Content") + theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 30),axis.text.y.left = element_text(size = 28),
                                                                            axis.title.y.left  = element_text(size = 30),axis.text.x = element_text(size = 28,angle=45,vjust = 0.4,hjust = 0.5),
                                                                            axis.ticks.length = unit(0.4,"cm"),axis.ticks = element_line(size = 1)) + xlab(label="")
PLOT4 <- ggplot(df.t[c(5091:10180),], aes(Method,Total, fill = Method)) + geom_violin(scale = "width") + geom_boxplot(width = 0.1,outlier.size = 0.5) + ggtitle(label="TIL10") +
  theme_classic() + ylab(label="Total Unnormalized Immune Content") + theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size = 30),axis.text.y.left = element_text(size = 28),
                                                                            axis.title.y = element_text(size = 30),axis.text.x = element_text(size = 28,angle=45,vjust = 0.4,hjust = 0.5),
                                                                            axis.ticks.length = unit(0.4,"cm"),axis.ticks = element_line(size = 1)) + xlab(label="")

#Save plot in PDF 10x36
ggsave(filename = "Figure_2/Figure_2e-h/Figure-_2e-f_Cell_Lines.pdf",plot=grid.arrange(PLOT1,PLOT2,PLOT3,PLOT4,ncol=4,nrow=1),device = "pdf",width = unit(40,"cm"),height = unit(10,"cm"))




