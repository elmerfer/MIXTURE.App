##This is a DEBUG VERSION


rm(list=ls())
library(data.table)
library(ComplexHeatmap)
library(ade4)#distancia de jaccard
library(ggplot2)
library(circlize)
library(dtangle)


## Change this path directory to the one where you download the code
##i.e .../.../my_directory/MIXTURE.R

.debug <- TRUE

##Load LM22 Signature
## Change this path directory to the one where you download the code
##i.e .../.../my_directory/LM22.RData
load("/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Current/LM22.RData")

##Utils
TestBetas <- function(mat, b.list){
  #This function estimate the root mean squared error between the estimated beta and the true simulated one
  #Args :
  # mat : the estimated betas (NxC) where N: number of samples, C: number of cell.types in the gene signature
  # b.list  : the list of simulated betas
  # 
  # Returns:
  # the RMSE for each mixture
  err <- vector("numeric", length(b.list))
  for ( i in 1: length(b.list)){
    err[i] <- sqrt(sum((mat[i,b.list[[i]]$id]-b.list[[i]]$beta[b.list[[i]]$id])^2))        
  }
  return(err)
}

DiffBetas <- function(mat, b.list){
  
  pred <- NULL
  bet <- NULL
  for ( i in 1: length(b.list)){
    pred <- c(pred,mat[i,b.list[[i]]$id] )
    bet <- c(bet,b.list[[i]]$beta[b.list[[i]]$id])
  }
  return(cbind(pred = pred, beta = bet))
}



TestExtraBetas <- function(mat, b.list){
  #This function estimate the extra proportion given by th extra coefficients
  #Args :
  # mat : the estimated betas (NxC) where N: number of samples, C: number of cell.types in the gene signature
  # b.list  : the list of simulated betas
  # 
  # Returns:
  # total extra proportion
  err <- vector("numeric", length(b.list))
  for ( i in 1: length(b.list)){
    err[i] <- sum(mat[ i, -b.list[[i]]$id] )        
  }
  return(err)
}


GetExtraBetas <- function(mat, b.list){
  bet <- NULL
  for( i in 1: length(b.list)){
    bet <- c(bet, mat[ i, -b.list[[i]]$id] )        
  }
  return(bet)
}

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
# write.table(M.pure.noise, file="/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/M.pure.mix.noise.xlsx", quote = FALSE, row.names = TRUE, sep="\t")
# 
# saveRDS(betas.noise.list, file = "/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/betas.noise.list.rds")


##RUNS SIMULATED SCENARIOS ----

##Scenario S1.1 - Prediction of Lm22)----
source('~/Dropbox/IDEAS/cibersort/MIXTURE/Utils/MIXTURE.DEBUG_V0.1.R')
load("/home/elmer/Dropbox/IDEAS/cibersort/MIXTURE/Data/LM22.RData")

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

#out.abis.ncrfe <- MIXTURE(expressionMatrix = M.c, signatureMatrix =  LM22, functionMixture =  nc.rfe.rlm, useCores = 1L)

##upload CIBERSORT proportions from CIBERSORT site for scaneario a
cib.site <- ReadCibersortWebResults(file="/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/LM22.Pure.CIBERSORT.Output_Job1.csv", type = "csv")

##Amount of estimated cell types
ciber<- apply(cib.site>0,1, sum) ##verificar si es 1 o 2
mixture <- apply(GetCellTypes(out.mixture)>0,1, sum)
ls.fit <- apply(GetCellTypes(out.ls)>0, 1, sum)
dt.fit <- apply(out.dt$estimates>0, 1, sum)
abis.fit <- apply(GetMixture(out.abis)>0,1,sum)

# df <- data.frame(id = c(1:length(ciber), 1:length(mixture), 1:length(ls.fit), 1:length(dt.fit)),N = c(ciber, mixture, ls.fit, dt.fit), 
#                  Method = c(rep("CIBERSORT",length(ciber)),rep("MIXTURE",length(mixture)), rep("ABBAS",length(ls.fit)), rep("DTANGLE", length(dt.fit))) )

summary(cbind(CIBERSORT=ciber,MIXTURE=mixture, ABBAS=ls.fit, DTANGLE=dt.fit, ABIS = abis.fit))
# CIBERSORT        MIXTURE      ABBAS          DTANGLE        ABIS      
# Min.   :1.000   Min.   :1   Min.   :6.000   Min.   :22   Min.   :10.00  
# 1st Qu.:2.000   1st Qu.:1   1st Qu.:6.000   1st Qu.:22   1st Qu.:11.00  
# Median :2.500   Median :1   Median :7.000   Median :22   Median :12.00  
# Mean   :3.227   Mean   :1   Mean   :7.091   Mean   :22   Mean   :11.91  
# 3rd Qu.:4.750   3rd Qu.:1   3rd Qu.:8.000   3rd Qu.:22   3rd Qu.:13.00  
# Max.   :8.000   Max.   :1   Max.   :9.000   Max.   :22   Max.   :13.00 



round(100*cumsum(table(ciber))/22,2)
# 1      2      3      4      5      6      7      8 
# 22.73  50.00  68.18  72.73  81.82  86.36  95.45 100.00 
round(100*cumsum(table(mixture))/22,2)
# 1 
# 100 
round(100*cumsum(table(ls.fit))/22,2)
# 6      7      8      9 
# 36.36  63.64  90.91 100.00 
round(100*cumsum(table(dt.fit))/22,2)
# 22
# 100
round(100*cumsum(table(abis.fit))/22,2)
# 10     11     12     13 
# 4.55  45.45  59.09 100.00 
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
M.cib.site <- cib.site
M.nu.robust <- GetMixture(out.mixture, "proportion")
M.abbas <- GetMixture(out.ls, "proportion")
M.dt <- out.dt$estimates
M.abis <- GetMixture(out.abis)

#autocorrelation of cell-types from the LM22 signature
CorrLM22 <- cor(M.c)


m <- max(cbind(M.cib.site, M.nu.robust,M.abbas, M.dt,M.abis),na.rm=TRUE)

col.map <- colorRamp2(c(0,m), c("blue","red"))
col.map.robust <- colorRamp2(c(0,m), c("grey","red"))
#renaming the cell types for graphical 
cell.types.names <- c("BN","BM","PC","CD8","CD4N","CD4Mr","CD4Ma","FH","Tr","TGD","NKr","NKa","M","M0","M1","M2","Dr","Da","Mr","Ma","E","N")
cell.types.names11 <- c("B","B","PC","CD8","CD4","CD4","CD4","CD4","CD4","TGD","NK","NK","Mo","Ma","Ma","Ma","D","D","Mt","Mt","Eo","N")

colnames(M.cib.site) <- rownames(M.cib.site) <- cell.types.names
dimnames(M.nu.robust) <- dimnames(M.cib.site)
dimnames(M.abbas) <- dimnames(M.cib.site)
dimnames(M.dt) <- dimnames(M.cib.site)
dimnames(CorrLM22) <- dimnames(M.cib.site)
dimnames(M.abis) <- dimnames(M.cib.site)
##preparing the matrices for visualization
##define null coefficients
M.cib.site[M.cib.site == 0] <- NA
M.abbas[M.abbas == 0] <- NA
M.dt[M.dt == 0] <- NA
M.nu.robust[M.nu.robust <= 0] <- NA
M.abis[M.abis == 0] <- NA


##Heatmaps for each model
setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/CIBERSORT.eps")##we can manage better
Heatmap(M.ciber.aux, cluster_rows = FALSE, show_row_names = FALSE, cluster_columns = FALSE, column_title = "CIBERSORT",name = "CIBERSORT",col = c("blue","red")) 
dev.off()

setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/CIBERSORTweb.eps")##we can manage better
Heatmap(M.cib.site, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "CIBERSORT",name = "CIBERSORT",col = c("blue","red")) 
dev.off()


setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/ABBAS.eps")##we can manage better
Heatmap(M.abbas, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "ABBAS",name = "ABBAS",col = c("blue","red")) 
dev.off()

setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/MIXTURE.eps")##we can manage better
Heatmap(M.nu.robust, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = c("red", "blue")) 
dev.off()

setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/DTANGLE.eps")##we can manage better
Heatmap(M.dt, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "DTANGLE",name = "DTANGLE",col = c("blue", "red")) 
dev.off()


setEPS()

M.abis2 <- M.abis
M.abis2[M.abis2 < 0] <- NA
diag(M.abis2) <-0 
diag(M.abis2)  <- max(M.abis2, na.rm=T )

postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/ABIS.eps")##we can manage better
Heatmap(M.abis2, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "ABIS",name = "ABIS",col = c("blue", "red")) 
dev.off()


summary(cbind(CIBERSORT=as.numeric(M.cib.site), MIXTURE = as.numeric(M.nu.robust), ABBAS= as.numeric(M.abbas), 
              DTANGLE = as.numeric(M.dt), ABIS=as.numeric(M.abis)))
# CIBERSORT         MIXTURE        ABBAS           DTANGLE               ABIS           
# Min.   :0.0001   Min.   :1     Min.   :0.0000   Min.   :0.0005039   Min.   :-0.2063857  
# 1st Qu.:0.0001   1st Qu.:1     1st Qu.:0.0006   1st Qu.:0.0024606   1st Qu.:-0.0027461  
# Median :0.0001   Median :1     Median :0.0022   Median :0.0052534   Median : 0.0002077  
# Mean   :0.3098   Mean   :1     Mean   :0.1410   Mean   :0.0454545   Mean   : 0.0454545  
# 3rd Qu.:0.9994   3rd Qu.:1     3rd Qu.:0.0072   3rd Qu.:0.0123500   3rd Qu.: 0.0066041  
# Max.   :0.9999   Max.   :1     Max.   :0.9996   Max.   :0.8120759   Max.   : 1.1447889  
# NA's   :413      NA's   :462   NA's   :328

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
setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/AbbasColinearity.eps")##we can manage better
Heatmap(M.abbas.aux2, cluster_rows = FALSE, show_row_names = TRUE, 
        cluster_columns = FALSE, column_title = "ABBAS",name = "ABBAS",
        col = colorRamp2(c(0, 0.04), c("blue",  "red"))) 
dev.off()

##multicolinearity for DTANGLE
M.dt.aux <- M.dt
diag(M.dt.aux) <- NA
summary(as.numeric(M.dt.aux))
diag(M.dt.aux) <- max(as.numeric(M.dt.aux), na.rm=T)

setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/DTANGLEColinearity.eps")##we can manage better
Heatmap(M.dt.aux, cluster_rows = FALSE, show_row_names = TRUE, 
        cluster_columns = FALSE, column_title = "DTANGLE",name = "DTANGLE",
        col = colorRamp2(c(0, 0.25), c("blue",  "red"))) 
dev.off()
M.abis2 <- M.abis
M.abis2[M.abis2 <= 0] <- NA

diag(M.abis2) <- 0
diag(M.abis2) <- max(M.abis2, na.rm=T)
dimnames(M.abis2) <- dimnames(CorrLM22)


#annot <- HeatmapAnnotation(text = anno_text(colnames(LM22), rot = 45, just = "left", offset = unit(2, "mm")))
## FULL COMPLETE FIGURE 1
hp <- Heatmap(CorrLM22, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, 
        cluster_columns = FALSE, column_title = "Cor LM22",name = "Cor LM22",col = c("blue", "red"),
        show_heatmap_legend = FALSE) +
Heatmap(M.abbas.aux2, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE,  
        column_title = "ABBAS",name = "ABBAS",
        col = colorRamp2(c(0, 0.04), c("blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(M.abis2, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE, 
          column_title = "ABIS",name = "ABIS",col = c("blue", "red"),show_heatmap_legend = FALSE) +
  Heatmap(M.dt.aux, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE,
            column_title = "DTANGLE",name = "DTANGLE",
          col = colorRamp2(c(0, 0.25), c("blue",  "red")),show_heatmap_legend = FALSE)+
Heatmap(M.cib.site, cluster_rows = FALSE, show_row_names = FALSE, cluster_columns = FALSE, column_title = "CIBERSORT",name = "CIBERSORT",
        col = c("blue","red"),show_heatmap_legend = FALSE, show_column_names = FALSE) +
  Heatmap(M.nu.robust, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",
          col = c("red", "blue"),show_heatmap_legend = FALSE, show_column_names = FALSE) 
setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/Colinearity.eps", paper = "a4", horizontal = FALSE, width = 800)##we can manage better
print(hp)
dev.off()



pdf("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/Colinearity.pdf", paper = "a4r", width = 12)##we can manage better
print(hp)
dev.off()

jpeg("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/Colinearity.jpg", width = 600)##we can manage better
print(hp)
dev.off()

Heatmap(CorrLM22, cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE, column_title = "Cor LM22",name = "Cor LM22",col = c("blue", "red"),top_annotation = annot, top_annotation_height = unit(5, "cm")) +
Heatmap(M.cib.site, cluster_rows = FALSE, show_row_names = FALSE, cluster_columns = FALSE, column_title = "CIBERSORT",name = "CIBERSORT",col = c("blue","red")) +
  Heatmap(M.nu.robust, cluster_rows = FALSE, show_row_names = TRUE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = c("red", "blue")) 




##------------------------------------------------------------------------------------------------######
## Scenario 1.S2)----
## evaluating the estimation of betas
## we randomly chose from 2 to 8 cell-types from LM22
## then, random proportions are choosen for each drawn from a uniform distritribution [0,1] 
## normalized by the sum of betas 
## then, a new mixture is build by LMM * betas


betas.list <- readRDS("/home/elmer/Dropbox/IDEAS/cibersort/MIXTURE/Data/betas.list.rds")

## assigning to M.c matrix the simulated cell-types mixtures
M.pure <- do.call(cbind, lapply(betas.list, function(x) x$A))


##estimating betas by CIBERSORT
out.betas.cib <- MIXTURE(expressionMatrix = M.pure, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE.old, useCores = 10L)
##estimating betas by MIXTURE

out.betas.robust <- MIXTURE(expressionMatrix = M.pure, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = 10L)

out.betas.abbas <- MIXTURE(expressionMatrix = M.pure, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = 10L)

out.betas.dt <- dtangle(Y = log2(t(M.pure)) , reference = log2(t(LM22)))

out.betas.abis <- MIXTURE(expressionMatrix = M.pure, signatureMatrix =  LM22, functionMixture =  rlm.abis, useCores = 1L)

# out.betas.ncrfe.abis <- MIXTURE(expressionMatrix = M.pure, signatureMatrix =  LM22, functionMixture =  nc.rfe.rlm, useCores = 1L)

##Getting the estimated proportions (beta hat)
M.r.cib <- ReadCibersortWebResults("/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/MixPure.CIBERSORT.Output_Job5.csv", type = "csv")

M.r.robust <- GetMixture(out.betas.robust, "prop")
M.r.abbas <- GetMixture(out.betas.abbas, "prop")
M.r.dt <- out.betas.dt$estimates
M.r.abis <- GetMixture(out.betas.abis)
# M.r.nc.abis <- GetMixture(out.betas.ncrfe.abis)

Beta <- do.call(rbind, lapply(betas.list, function(x) x$beta))


##Root mean squared error between estimated and true simulated Betas
df.BlandAltman <- data.frame(betahat = c( as.numeric(M.r.cib),M.r.robust,M.r.abbas,M.r.dt),
                             betasim = rep(as.numeric(Beta),4),
                             Method = rep(c("CIBERSORT","MIXTURE","ABBAS","DTANGLE"), each=prod(dim(M.r.cib))))

df.BlandAltman$difs <- df.BlandAltman$betahat - df.BlandAltman$betasim

df.summary <- ddply(df.BlandAltman, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T))

p <- ggplot(df.BlandAltman, aes(betasim, difs)) + 
  geom_point(na.rm=TRUE) + 
  geom_hline(data=df.summary,aes(yintercept=round(mean,3))) +
    geom_hline(data=df.summary,aes(yintercept=round(mean+2*sd,3))) + 
    geom_hline(data=df.summary,aes(yintercept=round(mean-2*sd,3))) +
   facet_wrap(~Method)

print(p)

ggplot(df.BlandAltman, aes(difs)) +   geom_density(aes(colour = Method)) 


 facet_wrap(~Method)


##Diff betas only account for the simulated ones.
cib.pred <- DiffBetas(M.r.cib, betas.list)
rob.pred <- DiffBetas(M.r.robust, betas.list)
abbas.pred <- DiffBetas(M.r.abbas, betas.list)
dt.pred <- DiffBetas(M.r.dt, betas.list)
summary(cbind(CIBERSORT=cib.pred[,1]-cib.pred[,2], MIXTURE= rob.pred[,1]-rob.pred[,2], 
              ABBAS= abbas.pred[,1]-abbas.pred[,2], DTANGLE=dt.pred[,1]-dt.pred[,2]))



par(mfrow=c(1,2))
boxplot(Difs)

summary(Difs)
t.test(Difs[,1], Difs[,2], paired = TRUE)

pairwise.wilcox.test(Difs, g = rep(c("CIBERSORT","MIXTURE","ABBAS"),each=1000), paired = TRUE)
# Pairwise comparisons using Wilcoxon signed rank test 
# 
# data:  Difs and rep(c("CIBERSORT", "MIXTURE", "ABBAS"), each = 1000) 
# 
# ABBAS  CIBERSORT
# CIBERSORT <2e-16 -        
#   MIXTURE   <2e-16 <2e-16   
# 
# P value adjustment method: holm

cibersort.web.length <- apply(M.r.cib>0, 1, sum)
robust.bet.length <- apply(GetCellTypes(out.betas.robust),1, sum)
abbas.bet.length <- apply(GetCellTypes(out.betas.abbas),1, sum)
dtangle.length <- apply(M.r.dt > 0, 1, sum)
abis.length <- apply(GetCellTypes(out.betas.abis),1, sum)
beta.length <- apply(Beta ,1, function(x) sum(x>0))

df.betas <- data.frame(est = c( abbas.bet.length, dtangle.length, cibersort.web.length, robust.bet.length,abis.length),
                       beta = as.factor(c(beta.length,beta.length,beta.length,beta.length,beta.length)),
                       method = factor(c(
                         rep("ABBAS",length(abbas.bet.length)),
                         rep("DTANGLE",length(dtangle.length)),
                         rep("CIBERSORT", length(cibersort.web.length)), 
                         rep("MIXTURE",length(robust.bet.length)),
                         rep("ABIS",length(abis.length))),
                         levels = c("ABBAS","ABIS","DTANGLE","CIBERSORT","MIXTURE")))


#Figure 1.C

setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/NumeberOfCoeffsPURE.eps")##we can manage better
p.betascount <- ggplot(df.betas, aes(x=beta, y=est, fill = method)) +
  geom_boxplot(position = position_dodge(1) ) + ylim(1,22) +
  labs( x = "True number of coefficients",
        y = "Estimated number of coefficients", 
        fill = "Algorithm")+ theme_bw() 
print(p.betascount)
dev.off()
PureNumBetasPlot <-   ggplot(df.betas, aes(x=beta, y=est, fill = method)) +
  geom_boxplot(position = position_dodge(1) ) + ylim(1,22) +
  labs( x = "True number of coefficients",
        y = "Estimated number of coefficients", 
        fill = "Algorithm") 
## error between expected betas and simulated ones
err.mixture <- TestBetas(GetMixture(out.betas.rfe, "prop"), betas.list)
err.robust <- TestBetas(GetMixture(out.betas.robust, "prop"), betas.list)
err.cibersort <- TestBetas(GetMixture(out.betas.cib, "prop"), betas.list)
err.abbas <- TestBetas(GetMixture(out.betas.cib, "prop"), betas.list)

df.betas.RMSE <- data.frame(est = c(err.cibersort,err.mixture, err.robust, err.abbas),
                       beta = as.factor(c(beta.length,beta.length,beta.length,beta.length)),
                       method = c(rep("CIBERSORT", length(cib.bet.length)), 
                                  rep("MIXTURE",length(rfe.bet.length)),
                                  rep("ROBUST",length(robust.bet.length)),
                                  rep("ABBAS",length(rfe.bet.length))))
df.betas.RMSE <- subset(df.betas.RMSE, method != "MIXTURE")


#Figure 1.C
ggplot(df.betas.RMSE, aes(x=beta, y=est, fill=method)) +
  geom_violin(position=position_dodge(1), trim =TRUE,
              draw_quantiles = c(0.5)) 



summary(cbind(CIBERSORT= TestExtraBetas(GetMixture(out.betas.cib, "prop"), betas.list),
   RFE= TestExtraBetas(GetMixture(out.betas.rfe, "prop"), betas.list),
   ROBUST=TestExtraBetas(GetMixture(out.betas.robust, "prop"), betas.list),
    ABBAS=TestExtraBetas(GetMixture(out.betas.abbas, "prop"), betas.list) ))
#         CIBERSORT             RFE              ROBUST        ABBAS          
# Min.   :4.219e-05   Min.   :-6.983e-05   Min.   :0   Min.   :7.330e-06  
# 1st Qu.:1.181e-04   1st Qu.: 2.041e-06   1st Qu.:0   1st Qu.:3.625e-03  
# Median :1.496e-04   Median : 2.369e-05   Median :0   Median :9.198e-03  
# Mean   :1.610e-04   Mean   : 3.566e-05   Mean   :0   Mean   :1.353e-02  
# 3rd Qu.:1.914e-04   3rd Qu.: 5.790e-05   3rd Qu.:0   3rd Qu.:2.017e-02  
# Max.   :5.223e-04   Max.   : 2.520e-04   Max.   :0   Max.   :7.705e-02 





## evaluating the estimation of betas + noise
## ## evaluating the estimation of betas
## we randomly chose from 2 to 8 cell-types from LM22,
## then, random proportions are choosen for each drawn from a uniform distritribution [0,1] 
## normalized by the sum of betas 
## random selection of 547 gene expression for the whole signature LM22 matrix
## then, a new mixture is build by LMM * betas + the random gene expression of the 547 genes

##Scenario 1.S3)----



betas.noise.list <- readRDS(file = "/home/elmer/Dropbox/IDEAS/cibersort/MIXTURE/Data/betas.noise.list.rds")
M.pure.noise <- do.call(cbind, lapply(betas.noise.list, function(x) x$A))



#M <- riaz[rownames(riaz) %in% rownames(LM22),]

out.betas.noise.cib <- MIXTURE(expressionMatrix = M.pure.noise, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE.old, useCores = 10L)
summary(as.vector(GetMixture(out.betas.noise.cib)))

out.robust.noise <- MIXTURE(expressionMatrix = M.pure.noise, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = 11L)

out.betas.noise.abbas <- MIXTURE(expressionMatrix = M.pure.noise, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = 3L)

dt_out.noise <- dtangle(Y=log2(t(M.pure.noise)) , reference = log2(t(LM22)))

out.betas.noise.abis <- MIXTURE(expressionMatrix = M.pure.noise, signatureMatrix =  LM22, functionMixture =  rlm.abis, useCores = 1L)

#out.betas.noise.ncrfe.abis <- MIXTURE(expressionMatrix = M.pure.noise, signatureMatrix =  LM22, functionMixture =  nc.rfe.rlm, useCores = 1L)


##from CIBERSORT site
M.r.cib.noise <-ReadCibersortWebResults("/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/NoiseMix.CIBERSORT.Output_Job4.csv", type = "csv")


Beta.noise <- do.call(rbind, lapply(betas.noise.list, function(x) x$beta))


cib.pred.n <- DiffBetas(M.r.cib.noise, betas.noise.list)
rob.pred.n <- DiffBetas(GetMixture(out.robust.noise, "prop"), betas.noise.list)
abbas.pred.n <- DiffBetas(GetMixture(out.betas.noise.abbas, "prop"), betas.noise.list)
dt.pred.n <- DiffBetas(dt_out.noise$estimates, betas.noise.list)
abis.pred.n <- DiffBetas(GetMixture(out.betas.noise.abis), betas.noise.list)



summary(cbind(CIBERSORT=cib.pred.n[,1]-cib.pred.n[,2], MIXTURE= rob.pred.n[,1]-rob.pred.n[,2], 
              ABBAS= abbas.pred.n[,1]-abbas.pred.n[,2], DTANGLE=dt.pred.n[,1]-dt.pred.n[,2], ABIS= abis.pred.n[,1]-abis.pred.n[,2]))
# CIBERSORT            MIXTURE               ABBAS              DTANGLE              ABIS          
# Min.   :-0.157742   Min.   :-0.1126468   Min.   :-0.442820   Min.   :-0.67744   Min.   :-0.202205  
# 1st Qu.:-0.016191   1st Qu.:-0.0051246   1st Qu.:-0.068592   1st Qu.:-0.15807   1st Qu.:-0.049209  
# Median :-0.008495   Median :-0.0015208   Median :-0.037011   Median :-0.10217   Median :-0.024475  
# Mean   :-0.011264   Mean   :-0.0027311   Mean   :-0.039437   Mean   :-0.11347   Mean   :-0.021731  
# 3rd Qu.:-0.003250   3rd Qu.: 0.0004756   3rd Qu.:-0.007921   3rd Qu.:-0.05695   3rd Qu.:-0.000216  
# Max.   : 0.167414   Max.   : 0.0516431   Max.   : 0.386868   Max.   : 0.84362   Max.   : 0.311339 

df.n.BlandAltman <- data.frame(betahat = c( as.numeric(M.r.cib.noise),
                                            as.numeric(GetMixture(out.robust.noise, "prop")),
                                            as.numeric(GetMixture(out.betas.noise.abbas, "prop")),
                                            as.numeric(dt_out.noise$estimates), as.numeric(GetMixture(out.betas.noise.abis))),
betasim = rep(as.numeric(Beta.noise),5),
                             Method = factor(rep(c("CIBERSORT","MIXTURE","ABBAS","DTANGLE","ABIS"), each=prod(dim(M.r.cib)))
                                             ,levels = c("ABBAS","ABIS","DTANGLE","CIBERSORT","MIXTURE")))

df.n.BlandAltman$difs <- df.n.BlandAltman$betahat - df.n.BlandAltman$betasim

df.n.summary <- ddply(df.n.BlandAltman, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
                      min = min(difs, na.rm=T), max= max(difs,na.rm=T))
df.n.summary
# Method          mean          sd        min       max
# 1     ABBAS  4.174317e-20 0.045549869 -0.4428203 0.4854852
# 2      ABIS -8.123358e-20 0.048755973 -0.2022052 0.4795753
# 3   DTANGLE -7.007989e-20 0.101957879 -0.6774381 0.9415084
# 4 CIBERSORT -2.063736e-08 0.011458019 -0.1577419 0.1674136
# 5   MIXTURE  9.543286e-19 0.005396148 -0.1126468 0.1069340
ddply(subset(df.n.BlandAltman, betasim ==0), .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
                      min = min(difs, na.rm=T), max= max(difs,na.rm=T))
# Method        mean          sd           min       max
# 1     ABBAS 0.011611194 0.032317013  0.0000000000 0.4854852
# 2      ABIS 0.006397955 0.047119913 -0.1638348817 0.4795753
# 3   DTANGLE 0.033408455 0.069108288  0.0001124162 0.9415084
# 4 CIBERSORT 0.003316475 0.007464299  0.0000000000 0.1528796
# 5   MIXTURE 0.000804109 0.004062058  0.0000000000 0.1069340
ddply(subset(df.n.BlandAltman, betasim >0), .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
      min = min(difs, na.rm=T), max= max(difs,na.rm=T))
# Method         mean          sd        min        max
# 1     ABBAS -0.039437222 0.059681059 -0.4428203 0.38686797
# 2      ABIS -0.021730545 0.047948859 -0.2022052 0.31133927
# 3   DTANGLE -0.113471243 0.113210436 -0.6774381 0.84362471
# 4 CIBERSORT -0.011264440 0.014958073 -0.1577419 0.16741360
# 5   MIXTURE -0.002731142 0.007894948 -0.1126468 0.05164312

p <- ggplot(df.n.BlandAltman, aes(betasim, difs)) + 
  geom_point(size =0.5) + 
  geom_hline(data=df.n.summary,aes(yintercept=round(mean,3)), color = "red", linetype = "solid") +
  geom_hline(data=df.n.summary,aes(yintercept=round(mean+2*sd,3)), color = "red", linetype = "dashed") + 
  geom_hline(data=df.n.summary,aes(yintercept=round(mean-2*sd,3)), color = "red", linetype = "dashed") +
  labs( x = "True Simulated coefficients",y = "Estimated - Simulated coefficients")+
  theme_bw()+ ggtitle("A")+ ylim(-.45, 0.6)+
  facet_wrap(~Method, nrow = 1)
setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/BlandAltmanEstimatedBetasNoise.eps")##we can manage better
print(p)
dev.off()

df.pn.summary <- ddply(df.n.BlandAltman[df.n.BlandAltman$betasim>0,], .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T))
ggplot(df.n.BlandAltman[df.n.BlandAltman$betasim>0,], aes(betasim, difs)) + 
  geom_point(size = 0.5) 
  geom_hline(data=df.pn.summary,aes(yintercept=round(mean,3))) +
  geom_hline(data=df.pn.summary,aes(yintercept=round(mean+2*sd,3))) + 
  geom_hline(data=df.pn.summary,aes(yintercept=round(mean-2*sd,3))) +
  facet_wrap(~Method)

sf <- subset(df.n.BlandAltman[df.n.BlandAltman$betasim>0,])#, Method %in% c("CIBERSORT","MIXTURE"))
sf.summary <- ddply(sf, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T))
setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/BlandAltmanEstimatedBetasPureNoise.eps")##we can manage better
ggplot(sf, aes(betasim, difs), colour = c("red","yellow")) + 
  geom_point(size = 0.5) + geom_smooth( method = "loess", span=0.5, fill = "green", color = "red", size = 0.5) + ylim(-0.13,0.08) +
  geom_hline(data=sf.summary,aes(yintercept=round(mean,3)), size =0.5) +
  geom_hline(data=sf.summary,aes(yintercept=round(mean+2*sd,3)), size =0.5) + 
  geom_hline(data=sf.summary,aes(yintercept=round(mean-2*sd,3)), size =0.5) +
  geom_hline(yintercept = 0, size =0.5, color = "blue") + 
  labs( x = "True Simulated coefficients",y = "Estimated - Simulated")+
  facet_wrap(~Method)
dev.off()


##bias for B>0
ddply(subset(df.n.BlandAltman, betasim > 0), .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), max = max(difs, na.rm=T), min = min(difs, na.rm=T))
# Method         mean          sd        max        min
# 1     ABBAS -0.039437222 0.059681059 0.38686797 -0.4428203
# 2      ABIS -0.021730545 0.047948859 0.31133927 -0.2022052
# 3   DTANGLE -0.113471243 0.113210436 0.84362471 -0.6774381
# 4 CIBERSORT -0.011264440 0.014958073 0.16741360 -0.1577419
# 5   MIXTURE -0.002731142 0.007894948 0.05164312 -0.1126468

##Bias for B==0
ddply(subset(df.n.BlandAltman, betasim == 0), .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), max = max(difs, na.rm=T), min = min(difs, na.rm=T))
# Method        mean          sd       max           min
# 1     ABBAS 0.011611194 0.032317013 0.4854852  0.0000000000
# 2      ABIS 0.006397955 0.047119913 0.4795753 -0.1638348817
# 3   DTANGLE 0.033408455 0.069108288 0.9415084  0.0001124162
# 4 CIBERSORT 0.003316475 0.007464299 0.1528796  0.0000000000
# 5   MIXTURE 0.000804109 0.004062058 0.1069340  0.0000000000
setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/EstimatedBetasNull.eps")##we can manage better
ggplot(subset(df.n.BlandAltman, betasim == 0), aes(Method, betahat, fill=Method)) + geom_boxplot()+
  labs( x = "Estimation Method",y = "Estimated coefficient for null simulated ones")
dev.off()


wilcox.test(difs~Method, subset(df.n.BlandAltman, betasim == 0 & Method %in% c("CIBERSORT","MIXTURE")), paired =TRUE)

# Wilcoxon signed rank test with continuity correction
# 
# data:  difs by Method
# V = 33319000, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
wilcox.test(difs~Method, subset(df.n.BlandAltman, betasim > 0 & Method %in% c("CIBERSORT","MIXTURE")), paired =TRUE)
# Wilcoxon signed rank test with continuity correction
# 
# data:  difs by Method
# V = 1004700, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

 
summary(cbind(CIBERSORT = cib.site.pred.n[,1]-cib.site.pred.n[,2], MIXTURE = rob.pred.n[,1]-rob.pred.n[,2], ABBAS = abbas.pred.n[,1]-abbas.pred.n[,2]))


##Comparison between simulated and estimated betas

Beta.noise<- do.call(rbind, lapply(betas.noise.list, function(x) x$beta))



beta.length.noise <- apply(Beta.noise ,1, function(x) sum(x>0))

##number of estimated betas
cib.bet.length.noise <- rowSums(M.r.cib.noise>0)

robust.bet.length.noise2 <- rowSums(GetMixture(out.robust.noise,"prop")>0)

abbas.bet.length.noise <- rowSums(GetMixture(out.betas.noise.abbas,"prop")>0)

dt.length.noise <- rowSums(dt_out.noise$estimates > 0) 

abis.bet.length.noise <- rowSums(GetMixture(out.betas.noise.abis,"prop")>0)

df.betas.noise <- data.frame(est = c(cib.bet.length.noise,
                                     robust.bet.length.noise2,
                                     abbas.bet.length.noise,
                                     dt.length.noise,
                                     abis.bet.length.noise),
                             beta = as.factor(c(beta.length.noise,beta.length.noise,beta.length.noise,beta.length.noise,beta.length.noise)),
                             method = factor(c(rep("CIBERSORT", length(cib.bet.length.noise)), 
                                      rep("MIXTURE",length(robust.bet.length.noise2)),
                                      rep("ABBAS",length(abbas.bet.length.noise)),
                                      rep("DTANGLE",length(beta.length.noise)),
                                      rep("ABIS",length(abis.bet.length.noise))),levels = c("ABBAS","ABIS","DTANGLE","CIBERSORT","MIXTURE")))

df.betas$Scenario <- "Pure"
df.betas.noise$Scenario <- "Noise"
df.total <- rbind(df.betas, df.betas.noise)
df.total$Scenario <- factor(df.total$Scenario, levels = c("Pure","Noise"))
##FIGURE 1
setEPS()
postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/FIG2.eps")##we can manage better
ggplot(df.total, aes(x=beta, y=est, fill = method)) +
  geom_boxplot(position = position_dodge(1) ) + ylim(1,22) +
  labs( x = "True number of coefficients",
        y = "Estimated number of coefficients", 
        fill = "Algorithm") + theme_bw() + facet_wrap(~Scenario) 
dev.off()


pdf("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/FIG2.pdf")##we can manage better
ggplot(df.total, aes(x=beta, y=est, fill = method)) +
  geom_boxplot(position = position_dodge(1) ) + ylim(1,22) +
  labs( x = "True number of coefficients",
        y = "Estimated number of coefficients", 
        fill = "Algorithm") + theme_bw() + facet_wrap(~Scenario) 
dev.off()

jpeg("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/FIG2.jpg")##we can manage better
ggplot(df.total, aes(x=beta, y=est, fill = method)) +
  geom_boxplot(position = position_dodge(1) ) + ylim(1,22) +
  labs( x = "True number of coefficients",
        y = "Estimated number of coefficients", 
        fill = "Algorithm") + theme_bw() + facet_wrap(~Scenario) 
dev.off()


