##This is a DEBUG VERSION
rm(list=ls())
library(data.table)
library(ComplexHeatmap)
library(ade4)#distancia de jaccard
library(ggplot2)
library(circlize)
library(dtangle)
library(immunedeconv)
source('Utils/MIXTURE.DEBUG_V0.1.R')

# source('~/Dropbox/IDEAS/cibersort/MIXTURE/Utils/my.quanTIseq.R')
## Change this path directory to the one where you download the code
##i.e .../.../my_directory/MIXTURE.R

.debug <- TRUE

##Load LM22 Signature
## Change this path directory to the one where you download the code
##i.e .../.../my_directory/LM22.RData
load("Data/LM22.RData")
load("Data/TIL10.RData")





##BUILD smulated scenarios a, b and c ----
##Scenario a)




# ##Scenario b)
# set.seed(123)
# betas.list <- lapply(1:1000, function(x, Mat, nrep) {
#   ns <-  sample(2:nrep,1)
#   r <- runif(ns, 0.2,1)
#   id <- sample(22,ns)
#   betas <- rep(0,22)
#   betas[id] <- r/sum(r)
#   A <- Mat %*% betas
#   list(beta = betas, id = id,A = A)
# }, data.matrix(LM22), nrep = 8)
# 
# ## assigning to M.c matrix the simulated cell-types mixtures
# M.pure <- do.call(cbind, lapply(betas.list, function(x) x$A))
# 
# 
# write.table(M.pure, file="/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/M.pure.mix.xlsx", quote = FALSE, row.names = TRUE, sep="\t")
# saveRDS(betas.list, file = "/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/betas.list.rds")
# 
# 
# ##Scenario c)
# set.seed(124)
# betas.noise.list <- lapply(1:1000, function(x, Mat, nrep) {
#   ns <-  sample(2:nrep,1)
#   r <- runif(ns, 0.2,1)
#   id <- sample(22,ns)
#   betas <- rep(0,22)
#   betas[id] <- r/sum(r)
#   A <- (Mat %*% betas) + matrix(as.vector(data.matrix(Mat))[sample(1:(nrow(Mat)*ncol(Mat)),nrow(Mat))],ncol=1) 
#   list(beta = betas, id = id,A = A)
# }, data.matrix(LM22), nrep = 8)
# 
# 
# 
# M.pure.noise <- do.call(cbind, lapply(betas.noise.list, function(x) x$A))
# write.table(M.pure.noise, file="Data/M.pure.mix.noise.xlsx", quote = FALSE, row.names = TRUE, sep="\t")
# 
# saveRDS(betas.noise.list, file = "Data/betas.noise.list.rds")
# set.seed(124)
# TIL10.betas.list <- lapply(1:1000, function(x, Mat, nrep) {
#   
#   ns <-  sample(2:nrep,1)
#   r <- runif(ns, 0.2,1)
#   id <- sample(ncol(Mat),ns)
#   betas <- rep(0,ncol(Mat))
#   betas[id] <- r/sum(r)
#   A <- (Mat %*% betas) 
#   list(beta = betas, id = id,A = A)
# }, Mat = data.matrix(TIL10), nrep = 6)
# 
# set.seed(124)
# TIL10.betas.noise.list <- lapply(1:1000, function(x, Mat, nrep) {
#   
#   ns <-  sample(2:nrep,1)
#   r <- runif(ns, 0.2,1)
#   id <- sample(ncol(Mat),ns)
#   betas <- rep(0,ncol(Mat))
#   betas[id] <- r/sum(r)
#   A <- (Mat %*% betas) + matrix(as.vector(data.matrix(Mat))[sample(1:(nrow(Mat)*ncol(Mat)),nrow(Mat))],ncol=1) 
#   list(beta = betas, id = id,A = A)
# }, Mat = data.matrix(TIL10), nrep = 6)

TIL10.M.pure <- do.call(cbind, lapply(TIL10.betas.list, function(x) x$A)) 

TIL10.M.pure.noise <- do.call(cbind, lapply(TIL10.betas.noise.list, function(x) x$A)) 

# saveRDS(TIL10.betas.list, file = "Data/TIL10.betas.list.rds")
# 
# saveRDS(TIL10.betas.noise.list, file = "Data/TIL10.betas.noise.list.rds")

# Save data to run CIBERSORT web
# aux <- data.frame("Gene Symbol" = rownames(TIL10.M.pure),TIL10.M.pure)
# write.table(aux, file="/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/TIL10.M.pure.txt", quote = FALSE, row.names = FALSE, sep="\t")
# aux <- data.frame("Gene Symbol" = rownames(TIL10.M.pure.noise),TIL10.M.pure.noise)
# write.table(aux, file="/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/TIL10.M.pure.noise.txt", quote = FALSE, row.names = FALSE, sep="\t")

##RUNS SIMULATED SCENARIOS ----

##Scenario S1.1 - Prediction of Lm22)----


M <- LM22
M.c <- M
##local implementation of cibersort 
#out.cib <- MIXTURE(expressionMatrix = M.c, signatureMatrix =  LM22, functionMixture =  cibersort, useCores = 3L)

##The CIBERSORT R version provides different results that its web site counterpart. so we dismiss it
#cib.dtang <- CIBERSORT(sig_matrix = LM22, mixture_file = M.c, perm = 0, QN = F) ## no produce los mismos resultados!!! from CIBERSORT_mod.R

out.mixture <- MIXTURE(expressionMatrix = M.c, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = 3L)

out.ls <- MIXTURE(expressionMatrix = M.c, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = 3L)

out.dt <- dtangle(Y = log2(t(LM22)) , reference = log2(t(LM22)))
##RLM
out.abis <- MIXTURE(expressionMatrix = M.c, signatureMatrix =  LM22, functionMixture =  rlm.abis, useCores = 1L)


out.quanti <-deconvolute_quantiseq.default(mix.mat = LM22, 
                                           arrays = FALSE, 
                                           signame = "/home/elmer/Dropbox/IDEAS/cibersort/MIXTURE/Data/LM22", 
                                           tumor = FALSE, 
                                           mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
# out.abis.ncrfe <- MIXTURE(expressionMatrix = M.c, signatureMatrix =  LM22, functionMixture =  nc.rfe.rlm, useCores = 1L)


##upload CIBERSORT proportions from CIBERSORT site for scaneario a
cib.site <- ReadCibersortWebResults(file="Data/LM22.Pure.CIBERSORT.Output_Job1.csv", type = "csv")

##TIL10
TIL10.out.mixture <- MIXTURE(expressionMatrix = TIL10, signatureMatrix =  TIL10, functionMixture =  nu.svm.robust.RFE, useCores = 3L)

TIL10.out.ls <- MIXTURE(expressionMatrix = TIL10, signatureMatrix =  TIL10, functionMixture =  ls.rfe.abbas, useCores = 3L)

##RLM
TIL10.out.abis <- MIXTURE(expressionMatrix = TIL10, signatureMatrix =  TIL10, functionMixture =  rlm.abis, useCores = 1L)


TIL10.out.quanti <-deconvolute_quantiseq.default(mix.mat = TIL10, 
                                                 arrays = FALSE, 
                                                 signame = "TIL10", 
                                                 tumor = FALSE, 
                                                 mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
TIL10.cib.web <-  ReadCibersortWebResults("Data/CIBERSORT.Output_TIL10.csv",type = "csv", nct = 10)

##Amount of estimated cell types LM22
ciber<- apply(cib.site>0,1, sum) ##verificar si es 1 o 2
mixture <- apply(GetCellTypes(out.mixture)>0,1, sum)
ls.fit <- apply(GetCellTypes(out.ls)>0, 1, sum)
abis.fit <- apply(GetMixture(out.abis)>0,1,sum)
quanti.fit <- apply(out.quanti[,-c(1,24)]>0,1,sum)

summary(cbind(CIBERSORT=ciber,MIXTURE=mixture, ABBAS=ls.fit, ABIS = abis.fit, QUANTISEQ = quanti.fit))
# CIBERSORT        MIXTURE      ABBAS            ABIS         QUANTISEQ    
# Min.   :1.000   Min.   :1   Min.   :6.000   Min.   :10.00   Min.   :1.000  
# 1st Qu.:2.000   1st Qu.:1   1st Qu.:6.000   1st Qu.:11.00   1st Qu.:1.000  
# Median :2.500   Median :1   Median :7.000   Median :12.00   Median :1.000  
# Mean   :3.227   Mean   :1   Mean   :7.091   Mean   :11.91   Mean   :1.182  
# 3rd Qu.:4.750   3rd Qu.:1   3rd Qu.:8.000   3rd Qu.:13.00   3rd Qu.:1.000  
# Max.   :8.000   Max.   :1   Max.   :9.000   Max.   :13.00   Max.   :2.000  

##amoount of estimated cell types for TIL10
TIL10.ciber<- apply(TIL10.cib.web>0,1, sum) ##verificar si es 1 o 2
TIL10.mixture <- apply(GetCellTypes(TIL10.out.mixture)>0,1, sum)
TIL10.ls.fit <- apply(GetCellTypes(TIL10.out.ls)>0, 1, sum)
TIL10.abis.fit <- apply(GetMixture(TIL10.out.abis)>0,1,sum)
TIL10.quanti.fit <- apply(TIL10.out.quanti[,-c(1,12)]>0,1,sum)
# df <- data.frame(id = c(1:length(ciber), 1:length(mixture), 1:length(ls.fit), 1:length(dt.fit)),N = c(ciber, mixture, ls.fit, dt.fit), 
#                  Method = c(rep("CIBERSORT",length(ciber)),rep("MIXTURE",length(mixture)), rep("ABBAS",length(ls.fit)), rep("DTANGLE", length(dt.fit))) )
summary(cbind(CIBERSORT=TIL10.ciber,MIXTURE=TIL10.mixture, ABBAS=TIL10.ls.fit, ABIS = TIL10.abis.fit, QUANTISEQ = TIL10.quanti.fit))
# CIBERSORT       MIXTURE      ABBAS           ABIS        QUANTISEQ
# Min.   :1.00   Min.   :1   Min.   :4.00   Min.   :4.00   Min.   :1  
# 1st Qu.:1.25   1st Qu.:1   1st Qu.:4.25   1st Qu.:4.00   1st Qu.:1  
# Median :2.00   Median :1   Median :5.00   Median :5.50   Median :1  
# Mean   :2.70   Mean   :1   Mean   :5.00   Mean   :5.80   Mean   :1  
# 3rd Qu.:4.00   3rd Qu.:1   3rd Qu.:5.75   3rd Qu.:7.75   3rd Qu.:1  
# Max.   :6.00   Max.   :1   Max.   :6.00   Max.   :8.00   Max.   :1 



##flating point error in CIBERSORT
cib.site.fpe <- cib.site
diag(cib.site.fpe) <- NA
cib.site.fpe[cib.site.fpe == 0] <- NA
summary(as.numeric(cib.site.fpe))

##fp errors for LS fit
ls.fpe <- GetMixture(out.ls,"prop")
diag(ls.fpe) <- NA
ls.fpe[ls.fpe == 0] <- NA
c(min(as.numeric(ls.fpe),na.rm=T),max(as.numeric(ls.fpe),na.rm=T))

##fpe DTANGLE
dt.fpe <- out.dt$estimates
diag(dt.fpe) <- NA
dt.fpe[dt.fpe == 0] <- NA
c(min(as.numeric(dt.fpe),na.rm=T),max(as.numeric(dt.fpe),na.rm=T))

##fpe in ABIS
df.abis <- GetMixture(out.abis)
diag(df.abis) <- NA
c(min(as.numeric(df.abis),na.rm=T),max(as.numeric(df.abis),na.rm=T))


##Required for FIGURE 1
##Estimated proportion for each estimated cell type i.e s_k, k!=i
##for LM22
M.cib.site <- cib.site
M.nu.robust <- GetMixture(out.mixture, "proportion")
M.abbas <- GetMixture(out.ls, "proportion")
M.dt <- out.dt$estimates
M.abis <- GetMixture(out.abis)
M.quanti <- data.matrix(out.quanti[,-c(1,24)])
#autocorrelation of cell-types from the LM22 signature
CorrLM22 <- cor(M.c)

#Estimated coefficientes for TIL10
TIL10.M.cib.site <- TIL10.cib.web
TIL10.M.nu.robust <- GetMixture(TIL10.out.mixture, "proportion")
TIL10.M.abbas <- GetMixture(TIL10.out.ls, "proportion")
TIL10.M.abis <- GetMixture(TIL10.out.abis)
TIL10.M.quanti <- data.matrix(TIL10.out.quanti[,-c(1,12)])
#autocorrelation of cell-types from the LM22 signature
CorrTIL10 <- cor(TIL10)


m <- max(cbind(M.cib.site, M.nu.robust,M.abbas, M.dt,M.abis),na.rm=TRUE)

col.map <- colorRamp2(c(0,m), c("blue","red"))
col.map.robust <- colorRamp2(c(0,m), c("grey","red"))
#renaming the cell types for graphical 
cell.types.names <- c("BN","BM","PC","CD8","CD4N","CD4Mr","CD4Ma","FH","Tr","TGD","NKr","NKa","M","M0","M1","M2","Dr","Da","Mr","Ma","E","N")
cell.types.names11 <- c("B","B","PC","CD8","CD4","CD4","CD4","CD4","CD4","TGD","NK","NK","Mo","Ma","Ma","Ma","D","D","Mt","Mt","Eo","N")

colnames(M.cib.site) <- rownames(M.cib.site) <- cell.types.names
dimnames(M.nu.robust) <- dimnames(M.cib.site)
dimnames(M.abbas) <- dimnames(M.cib.site)
# dimnames(M.dt) <- dimnames(M.cib.site)
dimnames(CorrLM22) <- dimnames(M.cib.site)
dimnames(M.abis) <- dimnames(M.cib.site)
dimnames(M.quanti) <- dimnames(M.cib.site)
##preparing the matrices for visualization
##define null coefficients
M.cib.site[M.cib.site == 0] <- NA
M.abbas[M.abbas == 0] <- NA
M.dt[M.dt == 0] <- NA
M.nu.robust[M.nu.robust <= 0] <- NA
M.abis[M.abis == 0] <- NA
M.quanti[M.quanti==0] <- NA

##For TIL10
TIL10.M.cib.site[TIL10.M.cib.site == 0] <- NA
TIL10.M.abbas[TIL10.M.abbas == 0] <- NA
TIL10.M.nu.robust[TIL10.M.nu.robust <= 0] <- NA
TIL10.M.abis[TIL10.M.abis == 0] <- NA
TIL10.M.quanti[TIL10.M.quanti==0] <- NA
TIL10.M.abis2 <- TIL10.M.abis
TIL10.M.abis2[TIL10.M.abis2 <= 0] <- NA

##FULL Figure TIL10 estimations
annot <- HeatmapAnnotation(text = anno_text(colnames(TIL10), rot = 45, just = "left", offset = unit(2, "mm")))
# FULL COMPLETE FIGURE 1
# hp <- Heatmap(CorrLM22, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE,
#         cluster_columns = FALSE, column_title = "Cor LM22",name = "Cor LM22",col = c("blue", "red"),
#         show_heatmap_legend = FALSE) +
TIL10.hp <- Heatmap(TIL10.M.abbas, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE,
                    column_title = "ABBAS",name = "ABBAS",
                    col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = T) +
  Heatmap(TIL10.M.abis2, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE,
          column_title = "ABIS",name = "ABIS",col = c("blue", "red"),show_heatmap_legend = FALSE) +
  Heatmap(TIL10.M.quanti, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE,
          column_title = "QUANTISEQ",name = "DTANGLE",col = c("blue", "red"),show_heatmap_legend = FALSE) +
  
  # Heatmap(M.dt.aux, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE,
  #           column_title = "DTANGLE",name = "DTANGLE",
  #         col = colorRamp2(c(0, 0.25), c("blue",  "red")),show_heatmap_legend = FALSE)+
  Heatmap(TIL10.M.cib.site, cluster_rows = FALSE, show_row_names = FALSE, cluster_columns = FALSE, column_title = "CIBERSORT",name = "CIBERSORT",
          col = c("blue","red"),show_heatmap_legend = FALSE, show_column_names = FALSE) +
  Heatmap(TIL10.M.nu.robust, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",
          col = c("red", "blue"),show_heatmap_legend = FALSE, show_column_names = FALSE)
setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/SelfTestTIL10.eps", paper = "a4", horizontal = FALSE, width = 800)##we can manage better
print(TIL10.hp)
dev.off()

##Heatmaps for each model
# setEPS()
# postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/CIBERSORT.eps")##we can manage better
# Heatmap(M.ciber.aux, cluster_rows = FALSE, show_row_names = FALSE, cluster_columns = FALSE, column_title = "CIBERSORT",name = "CIBERSORT",col = c("blue","red")) 
# dev.off()
# 
# setEPS()
# postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/CIBERSORTweb.eps")##we can manage better
# Heatmap(M.cib.site, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "CIBERSORT",name = "CIBERSORT",col = c("blue","red")) 
# dev.off()
# 
# 
# setEPS()
# postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/ABBAS.eps")##we can manage better
# Heatmap(M.abbas, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "ABBAS",name = "ABBAS",col = c("blue","red")) 
# dev.off()
# 
# setEPS()
# postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/MIXTURE.eps")##we can manage better
# Heatmap(M.nu.robust, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = c("red", "blue")) 
# dev.off()
# 
# setEPS()
# postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/DTANGLE.eps")##we can manage better
# Heatmap(M.dt, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "DTANGLE",name = "DTANGLE",col = c("blue", "red")) 
# dev.off()
# 
# setEPS()
#  postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/QUANTI.eps")##we can manage better
# Heatmap(M.quanti, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, 
#          column_title = "QUANTI",name = "DTANGLE",col = c("blue", "red")) 
#  dev.off()
# 
setEPS()

M.abis2 <- M.abis
M.abis2[M.abis2 < 0] <- NA
diag(M.abis2) <-0 
diag(M.abis2)  <- max(M.abis2, na.rm=T )

# postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/ABIS.eps")##we can manage better
# Heatmap(M.abis2, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "ABIS",name = "ABIS",col = c("blue", "red")) 
# dev.off()

summary(cbind(CIBERSORT=as.numeric(M.cib.site), MIXTURE = as.numeric(M.nu.robust), ABBAS= as.numeric(M.abbas), 
              DTANGLE = as.numeric(M.dt), ABIS=as.numeric(M.abis), QUANTI = as.numeric(M.quanti)))
# CIBERSORT         MIXTURE        ABBAS           DTANGLE               ABIS                QUANTI      
# Min.   :0.0001   Min.   :1     Min.   :0.0000   Min.   :0.0005039   Min.   :-0.2063857   Min.   :0.0044  
# 1st Qu.:0.0001   1st Qu.:1     1st Qu.:0.0006   1st Qu.:0.0024606   1st Qu.:-0.0027461   1st Qu.:0.8361  
# Median :0.0001   Median :1     Median :0.0022   Median :0.0052534   Median : 0.0002077   Median :1.0000  
# Mean   :0.3098   Mean   :1     Mean   :0.1410   Mean   :0.0454545   Mean   : 0.0454545   Mean   :0.8462  
# 3rd Qu.:0.9994   3rd Qu.:1     3rd Qu.:0.0072   3rd Qu.:0.0123500   3rd Qu.: 0.0066041   3rd Qu.:1.0000  
# Max.   :0.9999   Max.   :1     Max.   :0.9996   Max.   :0.8120759   Max.   : 1.1447889   Max.   :1.0000  
# NA's   :413      NA's   :462   NA's   :328                                               NA's   :458 

##range of fall detected celltypes


##The minimum ABBAS value is
min(as.numeric(M.abbas), na.rm=TRUE)
# [1] 2.01073e-05

#######################################################
## Preparing Figure 1 ########################

M.abbas.aux2 <- M.abbas
diag(M.abbas.aux2) <- 0.04
#the range of falsely detected  cell-types coefficient is [2e-5, 0.0323] 
#so we set the color limit to 0.04 to enhance this values in the heatmap.
#The color for the real coefficient will be "red" the full scale value of 0.04 spite
#the range was 0.0.9334 - 0.9669. It will be manually corrected in the Supplementary Fig.
# setEPS()
# postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/AbbasColinearity.eps")##we can manage better
# Heatmap(M.abbas.aux2, cluster_rows = FALSE, show_row_names = TRUE, 
#         cluster_columns = FALSE, column_title = "ABBAS",name = "ABBAS",
#         col = colorRamp2(c(0, 0.04), c("blue",  "red"))) 
# dev.off()

##multicolinearity for DTANGLE
M.dt.aux <- M.dt
diag(M.dt.aux) <- NA
summary(as.numeric(M.dt.aux))
diag(M.dt.aux) <- max(as.numeric(M.dt.aux), na.rm=T)

# setEPS()
# postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/DTANGLEColinearity.eps")##we can manage better
# Heatmap(M.dt.aux, cluster_rows = FALSE, show_row_names = TRUE, 
#         cluster_columns = FALSE, column_title = "DTANGLE",name = "DTANGLE",
#         col = colorRamp2(c(0, 0.25), c("blue",  "red"))) 
# dev.off()
M.abis2 <- M.abis
M.abis2[M.abis2 <= 0] <- NA

diag(M.abis2) <- 0
diag(M.abis2) <- max(M.abis2, na.rm=T)
dimnames(M.abis2) <- dimnames(CorrLM22)


# annot <- HeatmapAnnotation(text = anno_text(colnames(LM22), rot = 45, just = "left", offset = unit(2, "mm")))
## FULL COMPLETE FIGURE 1
# hp <- Heatmap(CorrLM22, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE,
# cluster_columns = FALSE, column_title = "Cor LM22",name = "Cor LM22",col = c("blue", "red"),
# show_heatmap_legend = FALSE) +
hp <- Heatmap(M.abbas.aux2, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE,
        column_title = "ABBAS",name = "ABBAS",
        col = colorRamp2(c(0, 0.04), c("blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(M.abis2, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE,
          column_title = "ABIS",name = "ABIS",col = c("blue", "red"),show_heatmap_legend = F) +
  Heatmap(M.quanti, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE,
            column_title = "QUANTISEQ",name = "col.scale",col = c("blue", "red"),show_heatmap_legend = T) +

  # Heatmap(M.dt.aux, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE,
  #           column_title = "DTANGLE",name = "DTANGLE",
  #         col = colorRamp2(c(0, 0.25), c("blue",  "red")),show_heatmap_legend = FALSE)+
Heatmap(M.cib.site, cluster_rows = FALSE, show_row_names = FALSE, cluster_columns = FALSE, column_title = "CIBERSORT",name = "CIBERSORT",
        col = c("blue","red"),show_heatmap_legend = FALSE, show_column_names = FALSE) +
  Heatmap(M.nu.robust, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",
          col = c("red", "blue"),show_heatmap_legend = F, show_column_names = FALSE)
#  setEPS()
#  postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/SelfTestLM22.eps", paper = "a4", horizontal = FALSE, width = 800)##we can manage better
  print(hp)
#  dev.off()
#  


  setEPS()
  postscript("BoxplotsEstimationNumber.eps")
  grid.arrange(hp, TIL10.hp, nrow = 2)
  dev.off()
  
# 
# 
# pdf("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/ColinearityLM22.pdf", paper = "a4r", width = 12)##we can manage better
# print(hp)
# dev.off()
# 
# jpeg("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/Colinearity.jpg", width = 600)##we can manage better
# print(hp)
# dev.off()
# 
# Heatmap(CorrLM22, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE, column_title = "Cor LM22",name = "Cor LM22",col = c("blue", "red"),top_annotation = annot, top_annotation_height = unit(5, "cm")) +
# Heatmap(M.cib.site, cluster_rows = FALSE, show_row_names = FALSE, cluster_columns = FALSE, column_title = "CIBERSORT",name = "CIBERSORT",col = c("blue","red")) +
#   Heatmap(M.nu.robust, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = c("red", "blue")) 

