library(openxlsx)
##Pls verify directories
source('Utils/MIXTURE.DEBUG_V0.1.R')
#Build data set for cibersort web 
##change your path and build your cellline.rds file. See Supplementary Material
cells <- readRDS("Data/celllines.rds")

dim(cells)

##get path to TIL10 signature file
load("Data/LM22.RData")
load("Data/TIL10.RData")
colnames(cells$targets)
M <- 2^cells$E #to run with MIXTURE and QUANTISEQ
cn[1:5]
cn <- colnames(cells$E)
rownames(M) <- cells$genes$symbol
# id <- c(0, seq(100,1000,100))
# for( i in 1:(length(id)-1) ){
#   M2 <- cbind(cells$genes$symbol, M[,(id[i]+1):id[i+1]])
#   colnames(M2) <- c("GeneSymbol", cn[(id[i]+1):id[i+1]])
#  write.table(M2, file = paste("/home/elmer/Elmer/MIXTURE/Cell",i,".txt",sep=""),sep="\t", quote=F, row.names = F)
# }
  #  M2 <- cbind(cells$genes$symbol, M[,1001:1018])
  #  colnames(M2) <- c("GeneSymbol", cn[1001:1018])
  # write.table(M2, file = paste("/home/elmer/Elmer/MIXTURE/Cell",11,".txt",sep=""),sep="\t", quote=F, row.names = F)

## RUN QUANTISEQ
# library(immunedeconv)
## We modify "deconvolute_quantiseq.default" in order to return both normalized and unnormalized coeficients
source('PaperMIXTURE/QUANTISEQ.test.R')
# 

#on LM22
#
# out.quanti <-deconvolute_quantiseq.default2(mix.mat = M, 
#                                             arrays = FALSE, 
#                                             signame = "/home/elmer/Dropbox/IDEAS/cibersort/MIXTURE/Data/LM22", 
#                                             tumor = FALSE, 
#                                             mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
 
 
cells.ls <- MIXTURE(expressionMatrix = M, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = 6L)
cells.ls.TIL10 <- MIXTURE(expressionMatrix = M, signatureMatrix =  TIL10, functionMixture =  ls.rfe.abbas, useCores = 6L)
# 
#  
#  ##RLM
  cells.abis <- MIXTURE(expressionMatrix = M, signatureMatrix =  LM22, functionMixture =  rlm.abis, useCores = 1L)
  cells.abis.TIL10 <- MIXTURE(expressionMatrix = M, signatureMatrix =  TIL10, functionMixture =  rlm.abis, useCores = 1L)

 # saveRDS(out.quanti, file=paste("/home/elmer/Elmer/MIXTURE/CellQUANTISEQ_LM22.RDS",sep=""))
#TIL10
# out.quanti <-deconvolute_quantiseq.default2(mix.mat = M, 
#                                             arrays = FALSE, 
#                                             tumor = FALSE, 
#                                             mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
# saveRDS(out.quanti, file=paste("/home/elmer/Elmer/MIXTURE/CellQUANTISEQ_TIL10.RDS",sep=""))

mix <- LoadMixtureResultsFromExcel("/home/elmer/Elmer/MIXTURE/outputCellines(2).xlsx")#Results from MIXTUREpy
MIX.Til10 <- LoadMixtureResultsFromExcel("/home/elmer/Elmer/MIXTURE/outputCellines_TIL10.xlsx")#Results from MIXTUREpy

#Load results from CIBERSORT
# cibA <- do.call(rbind, lapply(1:11, function(x){
# 
#   ci <- read.xlsx(paste("/home/elmer/Elmer/MIXTURE/CIBERSORT.Output_Abs_Cell",x,".xlsx",sep=""))
# }))
# #Load results from CIBERSORT
# cibA.Til10 <- do.call(rbind, lapply(1:11, function(x){
#   ci <- read.xlsx(paste("/home/elmer/Elmer/MIXTURE/CIBERSORT.Output_Abs_TIL10_Cell",x,".xlsx",sep=""))
# }))

# saveRDS(cibA,file=paste("/home/elmer/Elmer/MIXTURE/cibA.LM22.RDS",sep=""))
cibA <- readRDS(file=paste("/home/elmer/Elmer/MIXTURE/cibA.LM22.RDS",sep=""))
# saveRDS(cibA.Til10,file=paste("/home/elmer/Elmer/MIXTURE/cibA.Til10.RDS",sep=""))
cibA.Til10 <- readRDS(file=paste("/home/elmer/Elmer/MIXTURE/cibA.Til10.RDS",sep=""))
par(mfrow=c(1,2))
library(ggplot2)
#LM22
out.quanti <- readRDS(file=paste("/home/elmer/Elmer/MIXTURE/CellQUANTISEQ_LM22.RDS",sep=""))
colnames(cibA[,2:23])
df.lm22 <- data.frame(
  NumCels = c(rowSums(out.quanti$Abs[,-23]>0),
        rowSums(cibA[,2:23]>0),
        rowSums(GetMixture(mix, type = "abs")>0,na.rm=T),
        rowSums(GetMixture(cells.abis, type = "abs")>0,na.rm=T),
  rowSums(GetMixture(cells.ls, type = "abs")>0,na.rm=T)
  ),
  TotalI = c(rowSums(out.quanti$Abs[,-23]),
             rowSums(cibA[,2:23]),
             rowSums(GetMixture(mix, type = "abs"),na.rm=T),
             rowSums(GetMixture(cells.abis, type = "abs"),na.rm=T),
             rowSums(GetMixture(cells.ls, type = "abs"),na.rm=T)
  ),
  Method = factor(rep(c("QUANTISEQ","CIBERSORT","MIXTURE","ABIS","ABBAS"), each=nrow(GetMixture(mix))),
                  levels = c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")),
  TumorType = rep(cells$targets$cancer,5)
)





library(stringr)
ggplot(df.lm22, aes(Method, NumCels, fill = Method)) + geom_boxplot() + geom_jitter(,height = 0, width = 0.2)
ggplot(df.lm22, aes(Method, TotalI, fill = Method))  + geom_violin(scale = "width") + geom_boxplot(width = 0.1) + theme_bw()

ggplot(df.lm22[!str_detect(df.lm22$TumorType, "leukaemia"), ], aes(Method, NumCels, fill = Method)) + geom_boxplot() + geom_jitter(,height = 0, width = 0.2)
ggplot(df.lm22[!str_detect(df.lm22$TumorType, "leukaemia"), ], aes(Method, TotalI, fill = Method)) + geom_violin(scale = "width") + geom_boxplot(width = 0.1)#+ geom_jitter(width = 0.2, aes(shape = ".", colour ="grey") )


ggplot(subset(df.lm22,TumorType == "T_cell_leukemia"), aes(Method, TotalI, fill = Method)) + geom_violin() #+ geom_jitter(width = 0.2, aes(shape = ".", colour ="grey") )

ggplot(subset(df.lm22,Method == "MIXTURE" & TotalI > 0.04),  aes(TumorType, TotalI, fill = Method)) + 
  geom_violin() + geom_jitter(width = 0.2) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

kruskal.test(TotalI ~ Method, data = df.lm22)
wilcox.test(TotalI ~ Method, data = subset(df.lm22,Method!="CIBERSORT"), paired =TRUE)
wilcox.test(TotalI ~ Method, data = subset(df.lm22,Method!="MIXTURE"), paired =TRUE)
wilcox.test(TotalI ~ Method, data = subset(df.lm22,Method!="QUANTISEQ"), paired =TRUE)

wilcox.test(NumCels ~ Method, data = subset(df.lm22,Method!="CIBERSORT"), paired =TRUE)
wilcox.test(NumCels ~ Method, data = subset(df.lm22,Method!="MIXTURE"), paired =TRUE)
wilcox.test(NumCels ~ Method, data = subset(df.lm22,Method!="QUANTISEQ"), paired =TRUE)


summary(rowSums(cibA[,2:23]>0))
summary(rowSums(GetMixture(mix, type = "abs"),na.rm=T))
summary(rowSums(out.quanti$Abs[,-23]))



##TIL10
Q.Til10<-readRDS(file=paste("/home/elmer/Elmer/MIXTURE/CellQUANTISEQ_TIL10.RDS",sep=""))
colnames(Q.Til10$Abs)
df.til10 <- data.frame(
  NumCels = c(rowSums(Q.Til10$Abs[,-11]>0),
                rowSums(cibA.Til10[,2:11]>0),
                rowSums(GetMixture(MIX.Til10, type = "abs")>0,na.rm=T),
              rowSums(GetMixture(cells.abis.TIL10, type = "abs")>0,na.rm=T),
              rowSums(GetMixture(cells.ls.TIL10, type = "abs")>0,na.rm=T)
  ),
  TotalI = c(rowSums(Q.Til10$Abs[,-11]),
             rowSums(cibA.Til10[,2:11]),
             rowSums(GetMixture(MIX.Til10, type = "abs"),na.rm=T),
             rowSums(GetMixture(cells.abis.TIL10, type = "abs"),na.rm=T),
             rowSums(GetMixture(cells.ls.TIL10, type = "abs"),na.rm=T)
  ),
  Method = factor(rep(c("QUANTISEQ","CIBERSORT","MIXTURE","ABIS","ABBAS"), each=nrow(GetMixture(MIX.Til10))),levels = c("ABBAS","ABIS","QUANTISEQ","CIBERSORT","MIXTURE")),
  TumorType = rep(cells$targets$cancer,5)
)

ggplot(df.til10, aes(Method, NumCels, fill = Method)) + geom_boxplot() + geom_jitter(,height = 0, width = 0.2)
ggplot(df.til10, aes(Method, NumCels)) + geom_boxplot()# + geom_jitter()
ggplot(df.til10, aes(Method, TotalI, fill = Method)) + geom_violin(scale = "width") + geom_boxplot(width = 0.1) + theme_bw()

df.lm22$Signature <- "LM22"
df.til10$Signature <- "TIL10"
df.t <- rbind(df.lm22, df.til10)
colnames(df.t)

pdf("/home/elmer/Dropbox/IDEAS/cibersort/GenomeR/Figuras/CellsEstimates.pdf",paper="a4", width = 10, height = 3)
ggplot(df.t, aes(Method, TotalI, fill = Method)) + geom_violin(scale = "width") + geom_boxplot(width = 0.1,outlier.size = 0.5) + 
  theme_bw() + facet_grid(cols = vars(Signature))
dev.off()

pdf("/home/elmer/Dropbox/IDEAS/cibersort/GenomeR/Figuras/NumberCellsEstimates.pdf",paper="a4", width = 10, height = 3)
ggplot(df.t, aes(Method, NumCels, fill = Method)) +  geom_boxplot(,outlier.size = 0.0) + geom_jitter(height = 0, width = 0.2, size = 0.3) +
  theme_bw() + facet_grid(cols = vars(Signature))
dev.off()


ggplot(df.til10, aes(Method, NumCels)) + geom_violin() 
# ggplot(df.til10[!str_detect(df.til10$TumorType, "leukaemia"), ], aes(Method, TotalI)) +  geom_violin() 
# summary(df.til10[!str_detect(df.til10$TumorType, "leukaemia"), ]$Method)

kruskal.test(TotalI ~ Method, data = df.til10)
wilcox.test(TotalI ~ Method, data = subset(df.til10,Method!="CIBERSORT"), paired =TRUE)
wilcox.test(TotalI ~ Method, data = subset(df.til10,Method!="MIXTURE"), paired =TRUE)
wilcox.test(TotalI ~ Method, data = subset(df.til10,Method!="QUANTISEQ"), paired =TRUE)

wilcox.test(NumCels ~ Method, data = subset(df.til10,Method!="CIBERSORT"), paired =TRUE)
wilcox.test(NumCels ~ Method, data = subset(df.til10,Method!="MIXTURE"), paired =TRUE)
wilcox.test(NumCels ~ Method, data = subset(df.til10,Method!="QUANTISEQ"), paired =TRUE)

dim(MIX.Til10$Subjects$MIXabs)
dim(cibA.Til10)
boxplot(cbind(QUANTISEQ = rowSums(Q.Til10$Prop[,-c(1,12)]>0),
       CIBERSORT=rowSums(cibA.Til10[,2:11]>0),
      MIXTURE=rowSums(GetMixture(MIX.Til10, type = "pro")>0,na.rm=T)
))
apply(cbind(QUANTISEQ = rowSums(Q.Til10$Prop[,-c(1,12)]>0),
      CIBERSORT=rowSums(cibA.Til10[,2:11]>0),
      MIXTURE=rowSums(GetMixture(MIX.Til10, type = "pro")>0)
),2,summary)


apply(cbind(QUANTISEQ = rowSums(Q.Til10$Prop[,-c(1,12)]==0),
            CIBERSORT=rowSums(cibA.Til10[,2:11]==0),
            MIXTURE=rowSums(GetMixture(MIX.Til10, type = "pro")==0,na.rm=T)
),2,summary)
