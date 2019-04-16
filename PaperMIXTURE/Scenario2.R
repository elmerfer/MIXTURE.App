##REAL DATA from Noewman
library(dtangle)
source('Utils/MIXTURE.DEBUG_V0.1.R')
##since FL only contains B, CD8 and CD4 cells, we merge (sum) the proportions of such cell according to Newman et al.
cell.types.names11 <- c("B","B","PC","CD8","CD4","CD4","CD4","CD4","CD4","TGD","NK","NK","Mo","Ma","Ma","Ma","D","D","Mt","Mt","Eo","N")
CompositeMixture <- function(mixt, composite){
##This function merges and summarize the cell types proportion from the 22 cell-types to the 11 ones defined by Newman et al.  
  data.matrix(t(aggregate(data.frame(t(mixt)), by= list(composite), FUN = sum)[,-1]))
}


##load Data from Newman (downloaded from DTANGLE source code)
load("Data/newman_pbmc.rda")
load("Data/newman_fl.rda")

##prepare data for CIBERSORT and MIXTURE
FL <- t(2^newman_fl$data$log)[,1:14]
PBMC <- t(2^newman_pbmc$data$log)[,1:20]

##write data for CIBERSORT web site
# write.table(FL, file="/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/NewmanFL.txt", quote = FALSE, row.names = TRUE, sep="\t")
# write.table(PBMC, file="/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/NewmanPBCMC.txt", quote = FALSE, row.names = TRUE, sep="\t")

rn <- rownames(PBMC)

# for( i in 1:20){
#   df.out <- data.frame("Gene symbol" = rn)
#   df.out$pp <- PBMC[,i]
#   colnames(df.out) <- c("Gene symbol", paste("PBMC_",i,sep=""))
#   write.table(df.out, file=paste("/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/",colnames(df.out)[2],".csv",sep=""), quote = FALSE, row.names = FALSE, sep="\t")
# }


##Read Results from CIBERSORT web site
FL.cib.site <- ReadCibersortWebResults(file = "Data/FL_CIBERSORT.Output_Job6.csv")
PBMC.cib.site <- ReadCibersortWebResults(file = "Data/PBMC_CIBERSORT.Output_Job7.csv")

##read single subject CIBERSORT web

# FL.cib.ss <- do.call(rbind, lapply(1:14, function(x) {
#   print(paste("FL",x))
#   cat("\n")
#   read.csv(paste("/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/CIBres/FL",x,".csv", sep=""),h=T)
# } ))
# colnames(FL.cib.ss) <- str_replace_all(colnames(FL.cib.ss),"\\.", " " )
# colnames(FL.cib.ss)[10] <- "T cells regulatory  (Tregs)" 
# rownames(FL.cib.ss) <- paste("FL",1:14,sep=" ")
# 
# PBMC.cib.ss <- do.call(rbind, lapply(1:20, function(x) {
#   print(paste("PBMC",x))
#   read.csv(paste("/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/CIBres/PBMC",x,".csv", sep=""),h=T)
# } ))
# colnames(PBMC.cib.ss) <- str_replace_all(colnames(PBMC.cib.ss),"\\.", " " )
# colnames(PBMC.cib.ss)[10] <- "T cells regulatory  (Tregs)" 
# rownames(PBMC.cib.ss) <- paste("PBMC",1:20,sep=" ")


##Process by MIXTURE
FL.robust <- MIXTURE(expressionMatrix = FL, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = 5L)

PBMC.robust <- MIXTURE(expressionMatrix = PBMC, signatureMatrix =  LM22, functionMixture =  nu.svm.robust.RFE, useCores = 5L)

#LS fit
FL.abbas <- MIXTURE(expressionMatrix = FL, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = 3L)

PBMC.abbas <- MIXTURE(expressionMatrix = PBMC, signatureMatrix =  LM22, functionMixture =  ls.rfe.abbas, useCores = 5L)
#DTANGLE
dt_FL <- dtangle(Y=newman_fl$data$log[1:14,] , reference = newman_fl$data$log[-c(1:14),])

dt_PBMC <- dtangle(Y=newman_pbmc$data$log[1:20,] , reference = newman_pbmc$data$log[-c(1:20),])
##RLM
abis_FL <- MIXTURE(expressionMatrix = FL, signatureMatrix =  LM22, functionMixture =  rlm.abis, useCores = 1L)

abis_PBMC <- MIXTURE(expressionMatrix = PBMC, signatureMatrix =  LM22, functionMixture =  rlm.abis, useCores = 1L)


##single Subject

dt_FLss <- do.call(rbind,lapply(1:14,
                  function(x) dtangle(Y=newman_fl$data$log[x,, drop=F] , reference = newman_fl$data$log[-c(1:14),])$estimates
))




dt_PBMCss <- do.call(rbind,lapply(1:20,
                                function(x) dtangle(Y=newman_pbmc$data$log[x,, drop=F] , reference = newman_pbmc$data$log[-c(1:20),])$estimates
))

##Since FL data contains only B , CD8 and CD$ ct, we merge and summarise such CEll types according to Newman et al Supplementary File 
#https://www.nature.com/articles/nmeth.3337#supplementary-information
FL.cs.c <- CompositeMixture(FL.cib.site, cell.types.names11)
FL.cib.c.ss <-CompositeMixture(data.matrix(FL.cib.ss[,2:23]), cell.types.names11)
cbind(colnames(FL.cib.site),colnames(FL.cib.ss[,2:23]))
FL.r.c <- CompositeMixture(GetMixture(FL.robust,"prop"), cell.types.names11)
FL.a.c <- CompositeMixture(GetMixture(FL.abbas,"prop"), cell.types.names11)
FL.dt.c <- CompositeMixture(dt_FL$estimates, cell.types.names11)
FL.dt.c.ss <- CompositeMixture(dt_FLss, cell.types.names11)

FL.abis <- CompositeMixture(GetMixture(abis_FL,"prop"), cell.types.names11)
PBMC.abis <- CompositeMixture(GetMixture(abis_PBMC,"prop"), cell.types.names11)

class(FL.dt.c.ss)
sum(as.numeric(abs(FL.dt.c - FL.dt.c.ss)))

sum(as.numeric(abs(FL.cs.c - FL.cib.c.ss)))


dim(newman_fl$annotation$mixture[,c(3,4,2,1,5,6,7,8,9,10,11,12)])

df.FL <- data.frame(p = c(FL.cs.c[1,], FL.a.c[1,], FL.dt.c[1,], FL.r.c[1,],FL.abis[1,]) ,
                    model = factor( rep(c("CIBERSORT", "ABBAS", "DTANGLE","MIXTURE","ABIS"), 
      each = length(unique(cell.types.names11))), 
levels = c("ABBAS","ABIS","DTANGLE","CIBERSORT","MIXTURE")),
                  CT = rep(unique(cell.types.names11),times = 5),
                  truth = as.numeric(rep(newman_fl$annotation$mixture[1,c(3,4,2,1,5,6,7,8,9,10,11,12)], times = 5)), 
                  sj = rep(1, times = 12*5))



for( i in c(2:14)){
  df.FL <- rbind(df.FL, data.frame(p = c(FL.cs.c[i,], FL.a.c[i,], FL.dt.c[i,], FL.r.c[i,],FL.abis[i,]) ,
                                   model = factor( rep(c("CIBERSORT", "ABBAS", "DTANGLE","MIXTURE","ABIS"), 
                                                       each = length(unique(cell.types.names11))), 
                                                   levels = c("ABBAS","ABIS","DTANGLE","CIBERSORT","MIXTURE")),
                                   CT = rep(unique(cell.types.names11),times = 5),
                                   truth = as.numeric(rep(newman_fl$annotation$mixture[1,c(3,4,2,1,5,6,7,8,9,10,11,12)], times = 5)), 
                                   sj = rep(1, times = 12*5)))
}
df.FL$dif <- df.FL$p - df.FL$truth 


wilcox.test(dif~model, subset(df.FL, truth > 0 & model %in% c("CIBERSORT", "MIXTURE")), paired =TRUE, alternative = "less")

# Wilcoxon signed rank test with continuity correction
# 
# data:  dif by model
# V = 2408, p-value = 0.9998
# alternative hypothesis: true location shift is less than 0
wilcox.test(dif~model, subset(df.FL, truth > 0 & model %in% c("ABIS", "MIXTURE")), paired =TRUE, alternative = "less")
# Wilcoxon signed rank test
# 
# data:  dif by model
# V = 399, p-value = 0.2597
# alternative hypothesis: true location shift is less than 0
wilcox.test(dif~model, subset(df.FL, truth == 0 & model %in% c("CIBERSORT", "MIXTURE")), paired =TRUE, alternative = "greater")
# Wilcoxon signed rank test with continuity correction
# 
# data:  dif by model
# V = 2408, p-value = 0.0002182
# alternative hypothesis: true location shift is greater than 0
df.FL.noceros <- subset(df.FL, truth > 0)
ddply(subset(df.FL, truth == 0), .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T), max = max(dif, na.rm=T), min= min(dif, na.rm=T), q3 = quantile(dif, 0.95))
#true coefficientes == 0
# model       mean         sd       max          min        q3
# 1     ABBAS 0.04317701 0.07615312 0.2969441  0.000000000 0.2112136
# 2      ABIS 0.03280948 0.07327297 0.2410067 -0.103690014 0.1781785
# 3   DTANGLE 0.05247940 0.03966734 0.1660090  0.006456347 0.1253185
# 4 CIBERSORT 0.04067600 0.05614018 0.2110372  0.000000000 0.1651988
# 5   MIXTURE 0.03434851 0.06246182 0.2628134  0.000000000 0.1729970
# 
df.FL.summary <- ddply(df.FL, .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
df.FL.nc.summary <- ddply(df.FL.noceros, .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
#true coefficientes > 0
# model        mean        sd
# 1     ABBAS -0.12953102 0.1095302
# 2      ABIS -0.09842844 0.1124971
# 3   DTANGLE -0.15743819 0.2224409
# 4 CIBERSORT -0.12202703 0.1384909
# 5   MIXTURE -0.10304554 0.1223552
# 


p.FL <- ggplot(df.FL, aes(truth, dif)) + 
  geom_point(na.rm=TRUE, aes(colour = CT, shape = CT)) +
   geom_hline(data=df.FL.nc.summary,aes(yintercept=round(mean,3)), col = "red") +
   geom_hline(data=df.FL.nc.summary,aes(yintercept=round(mean+2*sd,3)), col = "red", linetype = "dashed") + 
   geom_hline(data=df.FL.nc.summary,aes(yintercept=round(mean-2*sd,3)), col = "red", linetype = "dashed") +
  geom_smooth(data = df.FL.noceros, span = 0.5) + 
  labs( x = "True Flow Cytometry Derived (FCD) proportions",y = "Error")  +
  theme_bw()+ theme(legend.position = "none") + ggtitle("B")+
  facet_wrap(~model, nrow = 1)


##Generate the figures#####
# setEPS()
# postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/NuevoFL_BlandAltman.eps")##we can manage better
# print(p.FL)
# dev.off()
# 
# png("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/NuevoFL_BlandAltman.png")##we can manage better
# print(p.FL)
# dev.off()



##PBMC Newman et al
##acomodar los nombres de las columnas para que sean iguales a los de LM22
library(stringr)
colnames(newman_pbmc$annotation$mixture) <- str_replace_all(colnames(newman_pbmc$annotation$mixture),"\\."," ")
colnames(newman_pbmc$annotation$mixture)
colnames(newman_pbmc$annotation$mixture)[which(str_detect(colnames(newman_pbmc$annotation$mixture), "T cells regulatory"))] <- "T cells regulatory (Tregs)"
newman_pbmc$annotation$mixture <- newman_pbmc$annotation$mixture[colnames(LM22)]

colnames(PBMC.cib.site)  <- str_replace_all( colnames(PBMC.cib.site), "\\." , " ")
colnames(PBMC.cib.site)[ which(str_detect(colnames(PBMC.cib.site), "T cells regulatory"))] <- "T cells regulatory (Tregs)"

xx <- as.numeric(data.matrix(newman_pbmc$annotation$mixture[1:20,]))
PBMC.cs.c <-PBMC.cib.site[,colnames(LM22)]
PBMC.r.c <- GetMixture(PBMC.robust,"prop")
PBMC.a.c <- GetMixture(PBMC.abbas,"prop")
PBMC.dt.c <- dt_PBMC$estimates[,colnames(LM22)]
PBMC.abis <- GetMixture(abis_PBMC,"prop")

df.t <- data.frame(x=c(xx,xx,xx,xx,xx), 
                   y = c(as.numeric(PBMC.cs.c),as.numeric(PBMC.r.c),as.numeric(PBMC.a.c),
                         as.numeric(PBMC.dt.c),as.numeric(PBMC.abis)),
                   model = factor(rep(c("CIBERSORT","MIXTURE","ABBAS","DTANGLE","ABIS"), each=length(as.numeric(PBMC.cs.c))),
                                  levels = c("ABBAS","ABIS","DTANGLE","CIBERSORT","MIXTURE")),
                   CT = rep(colnames(LM22),times = 5))
df.t$dif <- df.t$y-df.t$x
df.t.summary <- ddply(subset(df.t, x>0), .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
library(ggplot2)
df.PBMC <- df.t
df.PBMC.summary <- ddply(df.PBMC, .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
df.PBMC.summary <- ddply(subset(df.PBMC, x>0), .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
#true coefficientes > 0
# model         mean         sd
# 1     ABBAS -0.047056061 0.12019068
# 2      ABIS -0.019395463 0.15119402
# 3   DTANGLE -0.029325333 0.08241455
# 4 CIBERSORT -0.020569976 0.09303382
# 5   MIXTURE -0.009712779 0.09799701

ddply(subset(df.PBMC, x == 0), .(model), summarise, mean = mean(dif, na.rm=T), sd = sd(dif, na.rm=T))
#true coefficientes == 0
# model        mean         sd
# 1     ABBAS 0.032577273 0.06197797
# 2      ABIS 0.013427628 0.06609093
# 3   DTANGLE 0.020302154 0.02710329
# 4 CIBERSORT 0.014240776 0.02894889
# 5   MIXTURE 0.006724231 0.02589453

wilcox.test(dif~model, subset(df.PBMC, x > 0 & model %in% c("CIBERSORT", "MIXTURE")), paired =TRUE, alternative = "less")
# Wilcoxon signed rank test with continuity correction
# 
# data:  dif by model
# V = 4300, p-value = 0.007364
# alternative hypothesis: true location shift is less than 0
wilcox.test(dif~model, subset(df.PBMC, x > 0 & model %in% c("ABIS", "MIXTURE")), paired =TRUE, alternative = "less")
# Wilcoxon signed rank test with continuity correction
# 
# data:  dif by model
# V = 6529, p-value = 0.01051
# alternative hypothesis: true location shift is less than 0
wilcox.test(dif~model, subset(df.PBMC, x == 0 & model %in% c("CIBERSORT", "MIXTURE")), paired =TRUE, alternative = "greater")
# Wilcoxon signed rank test with continuity correction
# 
# data:  dif by model
# V = 5064, p-value = 4.368e-09
# alternative hypothesis: true location shift is greater than 0

col.cel.types <- c("chocolate1", "brown4", "black",
                   "tan1","green4", "green2", "lawngreen", "olivedrab3", "olivedrab", "chartreuse4",
                   "goldenrod","gold4","yellow","violetred","orangered1","red",
                   "plum4","plum","navy","mediumblue","cyan",
                   "grey28")
colores <- data.frame(CT=as.character(unique(df.PBMC$CT)), Colores=col.cel.types)
rownames(colores) <- colores[,1]

p.PBMC <-ggplot(df.PBMC, aes(x, dif)) + 
  geom_point(na.rm=TRUE, aes(colour = CT)) + 
  # scale_fill_manual(values  = as.character(colores[levels(df.t$CT),2]))+
  theme(legend.position = "none") +
  geom_hline(data=df.PBMC.summary,aes(yintercept=round(mean,3)),col="red") +
  geom_hline(data=df.PBMC.summary,aes(yintercept=round(mean+2*sd,3)), col = "red", linetype = "dashed") + 
  geom_hline(data=df.PBMC.summary,aes(yintercept=round(mean-2*sd,3)), col = "red", linetype = "dashed") + 
  geom_smooth(data = subset(df.t, df.t$x > 0), span =.7)+
  labs( x = "True Flow Cytometry Derived (FCD) proportions",y = "Error") +
  theme_bw()+ theme(legend.position = "none") + ggtitle("C")+
  facet_wrap(~model, nrow = 1)

print(p.PBMC)

# png("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/NewPBMC_BlandAltman.png")##we can manage better
# print(p.PBMC)
# dev.off()
# 
# setEPS()
# postscript("/home/elmer/Dropbox/IDEAS/cibersort/FiguresPaper/NewPBMC_BlandAltman.eps")##we can manage better
# print(p.PBMC)
# dev.off()

dens  <- lapply(c("CIBERSORT", "ABBAS", "DTANGLE","MIXTURE"), function(x) ggplot(subset(df.PBMC, model == x), aes(dif)) + geom_density() + coord_flip())
ba.plot <- lapply(c("CIBERSORT", "ABBAS", "DTANGLE","MIXTURE"), function(x) ggplot(subset(df.PBMC, model == x), aes(truth, dif)) + 
  geom_point(na.rm=TRUE) + 
  geom_hline(data=subset(df.PBMC.summary, model == "CIBERSORT"),aes(yintercept=round(mean,3))) +
  geom_hline(data=subset(df.PBMC.summary, model == "CIBERSORT"),aes(yintercept=round(mean+2*sd,3))) + 
  geom_hline(data=subset(df.PBMC.summary, model == "CIBERSORT"),aes(yintercept=round(mean-2*sd,3))))





