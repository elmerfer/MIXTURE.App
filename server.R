
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
# Authos: Elmer A. Fernandez
#

library(shiny)
# library(shinyjs)

library(DT)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(parallel)


source("Utils/MIXTURE.DEBUG_V0.1.R")

options(shiny.maxRequestSize=30*1024^2)

# .num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


# load("Data/LM22.RData")
# sigmat <- LM22

shinyServer( function(input, output, session) {
##INITIATLIZATION ####
  data <- reactiveValues( mixture.results = NULL , file2save = NULL)
  signature <- reactiveValues( Mat = NULL)
  dataset <- reactiveValues(sampleMat = NULL, filepath = NULL)
  
#   .Tabla <- NULL
  observeEvent(input$signature,{
    signature$Mat <- switch(input$signature,
                            LM22 = LM22,
                            TIL10 = TIL10)
  })
#   ######LOAD FUNCTIONS ###################

  observeEvent(input$mixResults, {
    if(str_detect(input$mixResults$datapath,"xlsx") == FALSE){
      showModal(modalDialog(
        title = "FILE ERROR",
        "Please choose an EXCEL file"
      ))
      return(NULL)
    }
    data$mixture.results <- LoadMixtureResultsFromExcel(input$mixResults$datapath)
    dataset$sampleMat <- NULL
    data$file2save <- NULL
    dataset$filepath <- NULL
    
    
      }  )

  ##called once you chose a file
  observeEvent(input$GeneExpr,{
    if(str_detect(input$GeneExpr$datapath,"xlsx") == FALSE){
      showModal(modalDialog(
        title = "FILE ERROR",
        "Please choose an EXCEL file"
      ))
      return(NULL)
    }else{
      dataset$filepath <- input$GeneExpr$datapath
      dataset$sampleMat <- read.xlsx(input$GeneExpr$datapath)
      rownames(dataset$sampleMat) <- dataset$sampleMat[,1]
      ##verify rownames 
      if( any(duplicated(dataset$sampleMat[,1]))) {
        ##hay duplicados
        m <- avereps(dataset$sampleMat[,-1] , ID = dataset$sampleMat[,1])
        dataset$sampleMat <-m
        # rownames(dataset$sampleMat) <- unique(dataset$sampleMat[,1])
            
      }else{
        
        dataset$sampleMat <- dataset$sampleMat[,-1]
        
      }
      
      data$file2save <- str_replace(dataset$name,".xlsx","MIXTURE_RES.xlsx")
      
      
    }
    
  })
  
  
  #   ###### CORE DECONVOLUTION #########
  ##I wnat to called them once i set the go button
  observeEvent(input$goButton, {
    # req(dataset())
    if( is.null(dataset$sampleMat)){
      showModal(modalDialog(
        title = "ERROR",
        "Please choose an Sample EXCEL file"
      ))
      return(NULL)  
    }
    
    

     # this would activate this render function
    if( is.null(signature$Mat) ){
      sigmat <- LM22
    }else{
      sigmat <- signature$Mat
    }
     progress <- shiny::Progress$new()
    
     
      nc <- input$cores
      it <- input$iter
      progress$set(message = "Computing: Wait until this window is closed", value=0)
      on.exit(progress$close())
      # if(input$MethodType == "M"){
        
        # genes <- intersect(rownames(dataset()), rownames(sigRNAseq))
        if(it >0){
         print("corriendo aca")
          data$mixture.results <- MIXTURE(expressionMatrix = dataset$sampleMat, signatureMatrix =  sigmat, functionMixture =  nu.svm.robust.RFE, useCores = nc,
                            nullDist = "PopulationBased", iter = it)
          
         
        }else{
          data$mixture.results <- MIXTURE(expressionMatrix = dataset$sampleMat, signatureMatrix =  sigmat, 
                                          functionMixture =  nu.svm.robust.RFE, useCores = nc) 
        }
        ##save results
      data$file2save <- str_replace(input$GeneExpr$name,".xlsx","_MIXTURE_RES.xlsx")
      SaveExcel(data$mixture.results, data$file2save)
  })
  #   ###### SUMMARY #########
    
    output$Summary <- renderPrint({
    nc <- input$cores
    it <- input$iter
    ngenes <- length(data$mixture.results$usedGenes)
    nsamples <- ncol(dataset$sampleMat)
    n.sample.genes <- nrow(dataset$sampleMat)
    if(is.null(dataset$sampleMat)){ 
      file.name <- "'No File selected'"
      }else{
        
        file.name <- input$GeneExpr$name
       ## output.file.name <- str_replace(input$GeneExpr$name,".xlsx", "_MIXTURE.xlsx")
        }
    if(is.null(signature$Mat)){
      sign.m <- "LM22"
    }else{
      # sign.m <- isolate(input$sigMat$name)
      sign.m <- isolate(input$signature)
    }
    cat("Number of cores: :",nc, " \nPermutations :" ,it, 
        "\nUsed Genes : ", ngenes, " / ",  n.sample.genes,
        "\nSignature Matrix :", sign.m,
        "\nFile :",file.name, " \nSamples :", nsamples)
    cat("\n")
    cat("Results file :", data$file2save)
  })
  #   ###### TABLES  #########
   output$table.proportions <- DT::renderDataTable({
     if(is.null(data$mixture.results)) return()
     df <- data$mixture.results$Subjects$MIXprop
     df[df==0] <- NA
     datatable(df) %>% formatRound(1:22, 3)
   },options = list(pageLength = 25, scrollX = TRUE))
  
   output$table.absolute <- DT::renderDataTable({
     if(is.null(data$mixture.results)) return()
     df <- data$mixture.results$Subjects$MIXabs
     df[df == 0 ] <- NA
     datatable(df) %>% formatRound(1:22, 3)
   },options = list(pageLength = 25))
   
   output$table.metrics <- DT::renderDataTable({
     if(is.null(data$mixture.results)) return()
     
     # df <- cbind(data$mixture.results$Subjects$ACCmetrix, pVal = data$mixture.results$p.values)
     df <- data$mixture.results$Subjects$ACCmetrix[,c(1,2,3,4,7,8,5,6)]
     datatable(df) %>% formatRound(1:6,3)
   },options = list(pageLength = 25))
   
   output$p.values <- DT::renderDataTable({
     if(is.null(data$mixture.results)) return()
     if(is.null(data$mixture.results$p.values)){
       showModal(modalDialog(
         title = "P values not available",
         "You should RUN MIXTURE with permutations"
       ))
       return()
     }
     # df <- cbind(data$mixture.results$Subjects$ACCmetrix, pVal = data$mixture.results$p.values)
     df <- data$mixture.results$p.values
  
     datatable(df) %>% formatRound(1:4,3)
   },options = list(pageLength = 25))
   
   output$used.genes <- DT::renderDataTable({
     if(is.null(data$mixture.results)) return()
     data.frame(Genes = data$mixture.results$usedGenes)
   },options = list(pageLength = 25))
   
  
   #   ###### SAVE #########
  output$downloadData <- downloadHandler(
      filename = function() {
        #data$file2save
          paste(str_remove(input$GeneExpr$name,".xlsx"),"_MIXTURE_.xlsx",sep="")
      },
      content = function(file) {
 
          SaveExcel(data$mixture.results, file)  
      
      },
      contentType = "Excel/xlsx"
    )
  
  
  
   #   ###### PLOTS #########
  output$bars <- renderPlot({
    if(is.null(data$mixture.results)) return()
    m.mix <- data$mixture.results$Subjects$MIXprop
    print(ncol(signature$Mat))
    df.test <- data.frame(b = as.vector(t(m.mix)), 
                          ct = rep(colnames(m.mix),nrow(m.mix)),
                          sbj = factor(rep(rownames(m.mix),each=ncol(signature$Mat)), levels = rownames(m.mix))) 
    col.cel.types <- c("chocolate1", "brown4", "black",
                       "tan1","green4", "green2", "lawngreen", "olivedrab3", "olivedrab", "chartreuse4",
                       "goldenrod","gold4","yellow","violetred","orangered1","red",
                       "plum4","plum","navy","mediumblue","cyan",
                       "grey28")
    col.cel.types <- col.cel.types[1:ncol(signature$Mat)]
    
    colores <- data.frame(CT=as.character(unique(df.test$ct)), Colores=col.cel.types)
    rownames(colores) <- colores[,1]
    
    ggplot(df.test, aes(sbj, b)) +   geom_col(aes(fill=ct)) + 
      scale_fill_manual(values  = as.character(colores[levels(df.test$ct),2])) + 
      theme(axis.text.x = element_text(angle = 90))+ 
      xlab("Subjects") + ylab("Proportions")
  })
  
  output$heatmap <- renderPlot({
    if(is.null(data$mixture.results)) return()
    df.ma <- data.frame(values = as.vector(t(data$mixture.results$Subjects$MIXprop)), 
                        CT = rep(colnames(data$mixture.results$Subjects$MIXprop),nrow(data$mixture.results$Subjects$MIXprop)),
                        Subjcets = rep(rownames(data$mixture.results$Subjects$MIXprop),each=ncol(data$mixture.results$Subjects$MIXprop)))
    df.ma$values[df.ma$values==0] <- NA
    ggplot(df.ma, aes(Subjcets, CT )) +
      geom_tile(aes(fill = values), color = "white") +
      scale_fill_gradient(low = "white", high = "red") +
      ylab("Cell Types ") +
      xlab("Subjects") +
      theme(legend.title = element_text(size = 10),
            legend.text = element_text(size = 12),
            plot.title = element_text(size=16),
            axis.title=element_text(size=14,face="bold"),
            axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(fill = "Proportion Levels")

  })
 
  output$immscoresuj <- renderPlot({
    if(is.null(data$mixture.results)) return()
    df <- data.frame(y = data$mixture.results$Subjects$ACCmetrix[,"IscBySbj"], 
                     x = rownames(data$mixture.results$Subjects$ACCmetrix))
    df$tipo = "In"
    mn <- mean(df$y, na.rm=T)
    mn.sd <- sd(df$y, na.rm=T)
    df$tipo[ df$y > mn + 2*mn.sd ] <- "Up"
    df$tipo[ df$y < mn + 2*mn.sd ] <- "Down"
    
    df.summary <- data.frame(mean = mean(df$y, na.rm=T), sd = sd(df$y, na.rm=T))
    
    ggplot(df, aes(x = x, y = y, colour = tipo)) + geom_point() + 
      geom_hline(data=df.summary,aes(yintercept=round(mean,3)), color = "red", linetype ="dashed") +
      geom_hline(data=df.summary,aes(yintercept=round(mean+2*sd,3)), color = "blue", linetype ="dashed") + 
      geom_hline(data=df.summary,aes(yintercept=round(mean-2*sd,3)), color = "blue", linetype ="dashed")+
      ylab("Immuno Score By Subject")+
      xlab("Subjects")+
      ggtitle("Subject based Immunmo Score factor") +
      theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1),legend.title = element_blank(), legend.position = "none")
    
  })
   
  output$immscorepob <- renderPlot({
    if(is.null(data$mixture.results)) return()
    df <- data.frame(y = data$mixture.results$Subjects$ACCmetrix[,"IscPob"], 
                     x = rownames(data$mixture.results$Subjects$ACCmetrix))
    df$tipo = "In"
    mn <- mean(df$y, na.rm=T)
    mn.sd <- sd(df$y, na.rm=T)
    df$tipo[ df$y > mn + 2*mn.sd ] <- "Up"
    df$tipo[ df$y < mn + 2*mn.sd ] <- "Down"
    
    df.summary <- data.frame(mean = mean(df$y, na.rm=T), sd = sd(df$y, na.rm=T))
    
    ggplot(df, aes(x = x, y = y, colour = tipo)) + geom_point() + 
      geom_hline(data=df.summary,aes(yintercept=round(mean,3)), color = "red", linetype ="dashed") +
      geom_hline(data=df.summary,aes(yintercept=round(mean+2*sd,3)), color = "blue", linetype ="dashed") + 
      geom_hline(data=df.summary,aes(yintercept=round(mean-2*sd,3)), color = "blue", linetype ="dashed")+
      ylab("Immuno score (Poblational)")+
      xlab("Subjects")+
      ggtitle("Population based Immunmo Score factor") +
      theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1),legend.title = element_blank(), legend.position = "none")
    
  })
  #   ###### ABOUT #########
  output$about <- renderPrint({
    cat("Author: : Elmer A. Fernández", " \nInstitution : CIDIE - UCC - CONICET (Argentina)" , "\nContact : efernandez at cidie.ucc.edu.ar")
  })
})
