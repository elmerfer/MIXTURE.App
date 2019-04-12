
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyjs)

library(parallel)
library(DT)
.num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
# sigmat <- LM22

navbarPage("MIXTURE",
           tabPanel("Run",
                    fluidRow(
                      column(4, wellPanel( h4("Files"),helpText("Signature Matrix is LM22 from Newman et al."),
                                           # fileInput('sigMat', strong('Choose Signature Matrix (optional):'), multiple = FALSE,
                                           #             accept = c('Excel', ".xlsx"), placeholder = "LM22 (Newman et al)"),
                                           
                                           
                                           fileInput('GeneExpr', strong('Choose gene expression file:'), multiple = FALSE,
                                                     accept = c('Excel', "xlsx") ),
                                           tags$hr(),
                                           fileInput("mixResults", strong('Load MIXTURE result file (EXCEL)'), multiple = FALSE,
                                                     accept = c('Excel', "xlsx"))
                      )),
                      column(4, wellPanel( h4("Options"),
                                           sliderInput("cores", "# CPUs", min = 1, max = .num.cores, value = .num.cores-1, step = 1 )
                       ),
                       wellPanel( h4("Save MIXTURE results (Excel)"),
                                    downloadButton("downloadData", "Download")
                       )  
                      ),
                      column(4, wellPanel( h4("Permutation Analysis"),
                                           # sliderInput("cores", "# CPUs", min = 1, max = .num.cores, value = .num.cores-1, step = 1 ),
                                           # Input: Specify the number of observations to view ----
                                           numericInput("iter", "Number of permutation samples:", value = 0,min = 0, step = 100),
                                           # 
                                          
                                           actionButton("goButton", "RUN MIXTURE")
                      ))
                      ),
                      fluidRow(
                      column(12, h4("Summary"),wellPanel(
                        verbatimTextOutput("Summary") 
                      ))  
                    )
                    ),
           navbarMenu("Tables",##tabPanel Tables
                      tabPanel("Proportions",
                               DT::dataTableOutput("table.proportions")),
                      tabPanel("Absolute",DT::dataTableOutput("table.absolute")),
                      tabPanel("Metrics",
                               fluidRow(
                                 column(12, wellPanel(DT::dataTableOutput("table.metrics")))#,
                                 #column(2)#,
                                 #column(4,  wellPanel(DT::dataTableOutput("used.genes")))
                                 
                               )
                               ),
                      tabPanel("P values",
                               fluidRow(
                                 column(12, wellPanel(DT::dataTableOutput("p.values")))#,
                                 #column(2)#,
                                 #column(4,  wellPanel(DT::dataTableOutput("used.genes")))
                               )
                      )
                      ),##fin menu Tables
           navbarMenu("Plots",##tabPanel Plots
                      tabPanel("Bars", plotOutput("bars")),
                      tabPanel("Heatmap", plotOutput("heatmap")),
                      tabPanel("Immuno Score By Subject", plotOutput("immscoresuj")),
                      tabPanel("Immuno Score Poulation Based", plotOutput("immscorepob"))
                      ),##Fin Menu Plots
           tabPanel("About",
                    fluidRow(
                      column(6, offset = 3,
                             verbatimTextOutput("about")
                      )))
)

