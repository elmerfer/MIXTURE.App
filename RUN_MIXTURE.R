##Script to RUN MIXTURE algorithm for LM22 deconvolution
##Autor: Elmer A. Fernandez (efernandez@cidie.ucc.edu.ar)
##Date: 28/11/2018
##-----------------------------------------------------------

## load MIXTURE functions
## Please verify the folder where you download the MIXTURE.R file
# source('~/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Current/MIXTURE/MIXTURE.R')
##Load signature matrix
library(openxlsx)
load("Data/LM22.RData")

##PLS VERIFY YOUR CURRENT DIRECTORY. IT SHOULD BE THE ONE WHERE YOU DOWNLOAD THE FILES.

source("Utils/MIXTURE.DEBUG_V0.1.R")

##Choose you sample file
sgm <- read.xlsx("Data/BRCA.subsample.xlsx")
## Verify if duplicated gene symbols
if( any(duplicated(sgm[,1]))){
  m <- avereps(sgm[,-1], ID=  sgm[,1])
  rownames(m ) <- unique(sgm[,1])
  sgm <- m
}else{
  rownames(sgm ) <- sgm[,1]
  sgm <- sgm[,-1]  
}
### multicore
## Verify your available cpu cores
num.cores <- 3L #if winfdows, only one is possible
##
mix.test <- MIXTURE(expressionMatrix = sgm,      #N x ncol(signatureMatrix) gene expresion matrix to evaluate 
                                                    ##rownames(M) should be the GeneSymbols
              signatureMatrix = LM22,               #the gene signature matrix (W) such that M = W*betas' (i.e the LM22 from Newman et al)
              iter = 0L,                            # amount of iteration in the statistical test (null distribution)
              functionMixture = nu.svm.robust.RFE,   #cibersort, nu.svm.optim.rfe, nnls = the cibersort model, 
                                                    ##nu-svm Recursive Feature Extraction and non negative lest squares
              useCores = num.cores,                        #cores for parallel processing
              verbose = TRUE,                       #TRUE or FALSE mesages  
              nullDist = "PopulationBased",                    #"none" or "PopulationBased", if the statistical test should be performed
              fileSave = "TETS_MIXTURE_RESULTS.xlsx")       #EXCEL file name to stare the results 

save(mix.test, file = "TEST_MIXTURE_FILE_LM22.RData") #save full lista as an RData object.
