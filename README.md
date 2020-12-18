![MIXTURE](https://github.com/elmerfer/MIXTURE.App/blob/master/www/Logo_B_1.pdf.png)

A v-SVR based noise constrained Recursive Feature Extraction algorithm for robust deconvolution of cell-types mixture from molecular signatures

Since the significant impact of immunotherapy in cancer, the estimation of the immune cell-type proportions present in a tumor becomes crucial. Currently, the deconvolution of the cell mixture content of a tumor is carried out by different analytic tools, yet the accuracy of inferred cell type proportions has room for improvement. We improve tumor immune environment characterization developing MIXTURE, an analytical method based on a noise constrained recursive variable selection for a support vector regression. 
[Get manuscript](https://www.biorxiv.org/content/10.1101/726562v1)

# NOTE: 

The MIXTURE shiny App has been only tested on Linux.
The RUN_MIXTURE code was tested on Linux, Windows and Mac. On windows only one CPU core is allowed.

# New! [MIXTURE in Python](https://github.com/MsMatias/MixturePy)
# New! [MIXTURE as an R library](https://github.com/elmerfer/MIXTURE)
# New! MIXTURE paper has been accepted in Brieffings in Bioinformatics! available soon!

## Test our Shiny app
You may test our shiny app on the free shiny server 
* You may use the following Excel [file](https://github.com/elmerfer/MIXTURE.App/blob/master/Data/BRCA.subsample.xlsx). Download it in your computer
* [Launch MIXTURE Shiny web app](https://cidie-conicet-ucc.shinyapps.io/mixture/)
* [How to use presentation](https://docs.google.com/presentation/d/1lv8YGpmyuf9n9UUKAm5GavVHrqdSYf9m1UrzU_a0sK8/edit?usp=sharing)


### Prerequisites

The current "functional like" version of the software requires the following libraries. (we recommend to use the [R library](https://github.com/elmerfer/MIXTURE) version)

* [data.table](https://cran.r-project.org/web/packages/data.table/)
* [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
* [ade4](https://cran.r-project.org/web/packages/ade4/index.html)
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
* [circlize](https://cran.r-project.org/web/packages/circlize/index.html)
* [e1071](https://cran.r-project.org/web/packages/e1071/index.html)
* [preprocessCore](https://www.bioconductor.org/packages/release/bioc/html/preprocessCore.html)
* [nnls](https://cran.r-project.org/web/packages/nnls/index.html)
* [plyr](https://cran.r-project.org/web/packages/plyr/index.html)
* [abind](https://cran.r-project.org/web/packages/abind/index.html)
* [openxlsx](https://cran.r-project.org/web/packages/openxlsx/index.html)
* [stringr](https://cran.r-project.org/web/packages/stringr/)


## Getting Started, Instllation and User Guide online instrucctions

Here you will find information regarding how to use the Shiny appplication as well as the comman line option. Please check it to see how to prepare your data file.

* [User Guide](https://docs.google.com/presentation/d/1lv8YGpmyuf9n9UUKAm5GavVHrqdSYf9m1UrzU_a0sK8/edit?usp=sharing)


## Running MIXTURE Shiny app, web Shiny app and command line

* [User Guide](https://docs.google.com/presentation/d/1lv8YGpmyuf9n9UUKAm5GavVHrqdSYf9m1UrzU_a0sK8/edit?usp=sharing)

This example tends to estimate the cell-types present in LM22 signature matrix from [Newman et al.](http://www.nature.com/nmeth/journal/v12/n5/abs/nmeth.3337.html) on some BRCA TCGA RNAseq samples
```
##PLS VERIFY YOUR CURRENT DIRECTORY. IT SHOULD BE THE ONE WHERE YOU DOWNLOAD THE FILES. 
## something like
##My favorite dir/Utils
##My favorite dir/Data

library(openxlsx)
load("Data/LM22.RData")



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

```

## Running on BRCA TCGA Data
Download the file BRCA_TCGA_MIXTURE_paper.R and open it. Follow the code. You will need to install the following packages:
* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
* [gplots](https://cran.r-project.org/web/packages/gplots/index.html)

The  TCGA BRCA data ready to use by the BRCA_TCGA_MIXTURE_paper.R script can be downloaded from [here](https://www.dropbox.com/s/zki1gkx5mq1quah/BRCA_rna.rds?dl=0)
## Contributing

Please read [CONTRIBUTING.md](https://github.com/elmerfer/MIXTURE.App/blob/master/Contributing.md) for details on our code of conduct, and the process for submitting pull requests to us.


## Authors

* **Elmer Andrés Fernández** - *Initial work* - [Profile](https://www.researchgate.net/profile/Elmer_Fernandez) - [CIDIE]- [CONICET](http://www.conicet.gov.ar) - [UCC](http://www.ucc.edu.ar)

## How to cite

Unveiling the immune infiltrate modulation in cancer and response to immunotherapy by MIXTURE—an enhanced deconvolution method
Elmer A Fernández, Yamil D Mahmoud, Florencia Veigas, Darío Rocha, Matías Miranda, Joaquín Merlo, Mónica Balzarini, Hugo D Lujan, Gabriel A Rabinovich, María Romina Girotti Briefings in Bioinformatics, bbaa317, https://doi.org/10.1093/bib/bbaa317 Published: 16 December 2020

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/elmerfer/MIXTURE.App/blob/master/LICENSE) file for details
