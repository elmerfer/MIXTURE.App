
library(data.table)
library(openxlsx)
library(org.Hs.eg.db)
library(limma)

#download and unzip expression from
url <- "https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/preprocessed/Cell_line_RMA_proc_basalExp.txt.zip"
fname <- basename(url)
download.file(url  =  url, destfile = fname, method = "auto")
unzip(fname)

#load expression matrix
a.data<- fread("Cell_line_RMA_proc_basalExp.txt")

#update annotation
b.annot<- a.data[,1:2]
colnames(b.annot)<- c("symbol", "name")
columns(org.Hs.eg.db)
b.entrezids <- mapIds(org.Hs.eg.db, keys=b.annot$symbol, column="ENTREZID", keytype="SYMBOL", multiVals="first")
b.entrezids[sapply(b.entrezids, is.null)]<- NA
b.annot$entrezid<- unlist(b.entrezids)

#fix colnames
colnames(a.data)<- gsub("DATA.", "", colnames(a.data))

#make elist
b.elist<- new("EList", list(E=a.data[,-c(1,2)], genes= b.annot))
dim(b.elist)

#remove missing entrezid
b.elist<- b.elist[which(!is.na(b.elist$genes$entrezid)),]

#combine repeated entrezid expression
b.elist<- avereps(x = b.elist, ID = b.elist$genes$entrezid)

saveRDS(b.elist, file = "Data/celllines.rds")
