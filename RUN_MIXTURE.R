##Script to RUN MIXTURE algorithm for LM22 deconvolution
##Autor: Elmer A. Fernandez (efernandez@cidie.ucc.edu.ar)
##Date: 28/11/2018
##-----------------------------------------------------------

## load MIXTURE functions
## Please verify the folder where you download the MIXTURE.R file
# source('~/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Current/MIXTURE/MIXTURE.R')
##Load signature matrix
load("~/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Current/LM22.RData")
sgm <- read.xlsx("/home/elmer/Dropbox/IDEAS/cibersort/MyCIBERTSORT/Debug/Melanoma78220.xlsx")
rownames(sgm ) <- sgm[,1]
sgm <- sgm[,-1]
sgm <- sgm[,1:5]
sgm <- read.xlsx("/home/elmer/Dropbox/IDEAS/cibersort/MIXTURE.App/BRCA.subsample.xlsx")
dim(sgm)
mix.test <- MIXTURE(expressionMatrix = sgm,      #N x ncol(signatureMatrix) gene expresion matrix to evaluate 
                                                    ##rownames(M) should be the GeneSymbols
              signatureMatrix = LM22,               #the gene signature matrix (W) such that M = W*betas' (i.e the LM22 from Newman et al)
              iter = 300,                            # amount of iteration in the statistical test (null distribution)
              functionMixture = nu.svm.robust.RFE,   #cibersort, nu.svm.optim.rfe, nnls = the cibersort model, 
                                                    ##nu-svm Recursive Feature Extraction and non negative lest squares
              useCores = 10L,                        #cores for parallel processing
              verbose = TRUE,                       #TRUE or FALSE mesages  
              nullDist = "PopulationBased",                    #"none" or "PopulationBased", if the statistical test should be performed
              fileSave = "TETS_MIXTURE_FILE_LM22.xlsx")       #EXCEL file name to stare the results 

save(mix.test, file = "MIXTURE_FILE_LM22.RData") #save full lista as an RData object.
