#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)
library(ggplot2)
library(shinythemes)
library(plotly)
library(DT)
library(shiny)
library(ggplot2)
library(shinythemes)
library(plotly)
library(shinycssloaders)
library(pbapply)
library(reshape)
library(reshape2)

# Define UI for application that draws a histogram
## ++1.1 read data in ################
readinBox <- fluidRow(titlePanel('Input GS data files'),
                          column(4,                        # structure of window
                             # 2st module
                             h4('GS genotypic data:'),
                             fileInput('markersData','Select the GS genotypic data file',
                                       accept = c(
                                           'text/csv',
                                           'text/comma-separated-values',
                                           'text/tab-separated-values',
                                           'text/plain',
                                           '.csv',
                                           '.tsv',
                                           '.txt'
                                       ),width = "70%"),
                             selectInput('markerFormat',"Input data format",c("Hapmap","Matrix","Custom"),width = "70%",selected = "Hapmap"),
                             checkboxInput("genome.headerCK", "Header", value = TRUE),
                             HTML('The GS markers accept hapmap format or matrix with sample in row and marker in column.')),
                      column(4,
                             ##phenotypic data
                             h4('GS phenotypic data:'),
                             fileInput('phenotypeData','Select the phenotypic data file',
                                       accept = c(
                                           'text/csv',
                                           'text/comma-separated-values',
                                           'text/tab-separated-values',
                                           'text/plain',
                                           '.csv',
                                           '.tsv',
                                           '.phe'
                                       ),width = "70%"),
                             selectInput('phenoFormat',"Input phenotype format",c("CSV","table text"),selected = "CSV",width = "70%"),
                             column(6,checkboxInput("phe.headerCK", "Header", value = TRUE)),
                             column(6,checkboxInput("phe.rownames", "Rownames", value = FALSE)),
                             HTML('The GS phenotype inputs accept csv or text with sample in row and trait in column.')),
                      
                      column(4,  
                             ##input subset of isoforms for investigation
                             h4('Other incidence matrices:'),
                             fileInput('file.Z','Select the fixed effects matrix data file',
                                       accept = c(
                                           'text/csv',
                                           'text/comma-separated-values',
                                           'text/tab-separated-values',
                                           'text/plain',
                                           '.csv',
                                           '.tsv'
                                       ),width = "70%"),
                             column(6,checkboxInput("fix.headerCK", "Header", value = FALSE)),
                             column(6,checkboxInput("fix.rownamesCK", "Row names", value = FALSE)),
                             HTML('The fixed effects incidence matrices for model.'))
                     
)

## ++1.2 data control ###########
dataQC  <- fluidRow(titlePanel('Data preprocessing'),
                     sidebarLayout(   
                         sidebarPanel(
                             h4(strong('Parameter setting')),
                             br(),
                             numericInput("filter_maf","Minor allele frequency ( < will be remove)",0.05,0,1,0.01),
                             hr(style = "border-top: 1px solid lightgrey;"),
                             numericInput("filte_mr","Missing rate ( > will be remove)",0.20,0,1,0.01),
                             hr(style = "border-top: 1px solid lightgrey;"),
                             # HTML('<b style="color:grey;">Press GSDataQC button to set the parameters.</b>'),      
                             # HTML('<br><br>'),
                             h4(strong('Imputation')),
                             selectInput('filter_imputation',"Do imputation?",c("FALSE","TRUE")),
                             conditionalPanel(
                               condition = "input.filter_imputation == 'TRUE'",
                               radioButtons('filter_imputation_method', label = 'Imputation method',c("Major" = "median","Mean" = "mean"),selected = "median",inline = T),
                              ),
                             actionButton('dataPreprocessing','Marker data processing',icon("send outline icon"),class="btn btn-primary",width = "100%")
                             ,
                             HTML('<b style="color:red;">Press above button for marker data pre-processing.</b>'),
                             textOutput("runProcessIdx"),
                             hr(style = "border-top: 1px solid lightgrey;"),
                             h4(strong('Download full data')),
                            
                             column(6,downloadButton('download.MAF.Info', 'MAF information',class="btn btn-primary")),
                             column(6,downloadButton('download.numeric.data', 'Marker matrix',class="btn btn-primary")),
                             
                             HTML('<b style="color:red;">Press above button to download full data.</b>'),
                             # style = "height:600px",
                             width = 3,
                             HTML("<p></p>")
                         ),
                         mainPanel(
                             # use_waiter(),
                             # waiter_show_on_load(spin_fading_circles()),
                             tabsetPanel(type = "tabs",
                                         tabPanel("Markers data",
                                                  HTML('<b>Table:</b> The raw data of GS markers, default shows the first 100 lines and 100 columns.<br/>'),
                                                  br(),
                                                  DTOutput("markerView") %>% withSpinner(color="#0dc5c1"),
                                                  style = "height:550px; overflow-y: scroll;overflow-x: scroll;"
                                                  
                                         ),
                                         tabPanel("Statistic for raw data",
                                                  HTML('<b>Table:</b> Statistic of genotype data.<br/>'),
                                                  br(),
                                                  DTOutput('sta_maf_Info') %>% withSpinner(color="#0dc5c1"),
                                                  ## define box options 
                                                  style = "height:550px; overflow-y: scroll;overflow-x: scroll;" 
                                         ),
                                         tabPanel("Processed marker data",
                                                  HTML('<b>Table:</b> The genotype data after data processing and transformation, default shows the first 100 lines and 100 columns.<br/>'),
                                                  br(),
                                                  DTOutput('pro_marker') %>% withSpinner(color="#0dc5c1"),
                                                  ## define box options 
                                                  style = "height:550px; overflow-y: scroll;overflow-x: scroll;" 
                                         ),
                                         tabPanel("Phenotypic data",
                                                  HTML('<b>Table:</b> The input data of phenotypes.<br/>'),
                                                  br(),
                                                  DTOutput('pheRaw') %>% withSpinner(color="#0dc5c1"),
                                                  ## define box options 
                                                  style = "height:550px; overflow-y: scroll;overflow-x: scroll;" 
                                         ),
                                         tabPanel("PCA data",
                                                  HTML('<b>Table:</b> The result data of PCA<br/>'),
                                                  br(),
                                                  DTOutput('PCARes') %>% withSpinner(color="#0dc5c1"),
                                                  ## define box options 
                                                  style = "height:550px; overflow-y: scroll;overflow-x: scroll;" 
                                         ),
                                         tabPanel("UAMP data",
                                                  HTML('<b>Table:</b> The result data of UAMP<br/>'),
                                                  br(),
                                                  DTOutput('UMAPRes') %>% withSpinner(color="#0dc5c1"),
                                                  ## define box options 
                                                  style = "height:550px; overflow-y: scroll;overflow-x: scroll;" 
                                         ),
                                         tabPanel("TSNE data",
                                                  HTML('<b>Table:</b> The result data of T-SNE <br/>'),
                                                  br(),
                                                  DTOutput('TSNERes') %>% withSpinner(color="#0dc5c1"),
                                                  ## define box options 
                                                  style = "height:550px; overflow-y: scroll;overflow-x: scroll;" 
                                         )
                             ),width = 9
                         )
                     )
)


## ++1.3 input data visulization ##########
inputDataVisual <- fluidRow(
  titlePanel('Visualization of input data'),
  sidebarLayout(
    sidebarPanel(
             h4(strong('Plot panel size setting')),
             sliderInput('plotPanelWidth', label = 'Width: ', min = 50, max = 100,value = 100,step = 1),
             sliderInput('plotPanelHeight', label = 'Height: ', min = 400, max = 1500,value = 800,step = 50),
      hr(style = "border-top: 1px solid lightgrey;"),
      fluidRow(
        column(12,
               radioButtons("plotInputData", label = h4("Select plot data:"),
                            choices = list("Genotype" = 'genotype',"Phenotype" = 'phenotype'),
                            selected = 'phenotype',inline = T)
        ),

        column(12,
        conditionalPanel(
          condition = "input.plotInputData == 'genotype'",
          # radioButtons("genotypePoltType", label = h4("Select plot type:"),
          #              choices = list("PCA" = 'PCA',"UMAP" = 'UMAP',"T-SNE" = "tsne"),
          #              selected = 'PCA',inline = T),
          checkboxGroupInput("genotypeComputType", label = h4("Select plot type:"),
                             choices = list("PCA" = 'PCA',"UMAP" = 'UMAP',"T-SNE" = "TSNE"),
                             selected = 'PCA',inline = T),
          numericInput('cutTreeCount', label = 'Group Numbers default: ',value  = 3),
          numericInput('dimensionSelect', label = 'Select first N dimension: ',value  = 3),
          actionButton("genotypeComputStart", label = "Start to structure the population",icon("send outline icon"),class="btn btn-primary",width = "100%"),
          HTML("<p></p>"),
          br(),
          # downloadButton('downloadComputRes', 'PCA/UMAP/T-SNE information ...',class="btn btn-primary",width = "200px"),
          uiOutput("genotypePlotType"),
          fileInput('groupInformation','Select the file of group labels and colors.',
                    accept = c(
                      'text/csv',
                      'text/comma-separated-values',
                      'text/tab-separated-values',
                      'text/plain',
                      '.csv',
                      '.tsv',
                      '.txt'
                    ),width = "100%"),
          # sliderInput("groupNumbers", "Group numbers:", min = 1, max = 10,value = 3,step = 1),
          textInput('genotypePointColor', label = 'Points color: ',value  = "red,blue,green"),
          h6("Default color palette of rainbow"),
          # sliderInput("genotypePointSize", "Points size:", min = 0, max = 10,value = 2,step = 0.1)
          # sliderInput('D3.color.high', label = 'Color high: ',min = 0,max = 8000,value = 7000)
        ),
        conditionalPanel(
          condition = "input.plotInputData == 'phenotype'",
          radioButtons("phenotypePlotType", label = h4("Select plot type:"),
                       choices = list("Histogram" = 'histogram',"Density" = 'density',"Pairs" = "pairs"),
                       selected = 'density',inline = T),
          conditionalPanel(
            condition = "input.phenotypePlotType == 'histogram' | input.phenotypePlotType == 'density'",
            uiOutput("traitSelectionPanel"),
            conditionalPanel(
              condition = "input.phenotypePlotType == 'histogram'",
              sliderInput("binsN", "Number of bins to view:", min = 0, max = 100,value = 10),
              textInput('binsColor', label = 'Box color: ',value  = "blue"),
              radioButtons("showDensityLine", label = h4("Density/frequency line:"),
                           choices = list("FALSE" = 'FALSE',"TRUE" = 'TRUE'),
                           selected = 'FALSE',inline = T),
              conditionalPanel(
                condition = "input.showDensityLine == 'TRUE'",
                textInput('densityLineColor', label = 'Density line color: ',value  = "red"),
                sliderInput('densityLineSize', label = 'Density line size: ', min = 0, max = 10,value = 2,step = 0.1)
              )
            ),
            conditionalPanel(
              condition = "input.phenotypePlotType == 'density'",
              textInput('densityLineColorD', label = 'Density line color: ',value  = "blue"),
              sliderInput('densityLineSizeD', label = 'Density line size: ', min = 0, max = 10,value = 2,step = 0.1)
            )
          ),
        )
        )
        
      ), width = 3
    ),
     mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Visulization of genotype",
                           # HTML('<b>Table:</b> The raw data of phenotypic data'),
                           # plotlyOutput("genotypePlot",height = "800px"),
                           uiOutput("genotypePlot") %>% withSpinner(color="#0dc5c1")
                           # style = "height:580px; overflow-y: scroll;overflow-x: scroll;"
                  ),
                  tabPanel("Visulization of phenotype",
                           # HTML('<b>Table:</b> Statistic of genotype data'),
                           # br(),
                           # plotlyOutput('phenotypePlot',height = "800px"),
                           uiOutput("phenotypePlot") %>% withSpinner(color="#0dc5c1")
                           ## define box options
                           # style = "height:580px; overflow-y: scroll;overflow-x: scroll;"
                  )
                 ), width = 9
      )
    )
    
)


## 1. page 1 #############
page1 <- tabPanel("G2P data input and quality control", #tabPanel the 1st module and title
                  fluidPage(
                      readinBox,##input data
                      hr(style = "border-top: 1px solid lightgrey;"),
                      
                      dataQC,
                      hr(style = "border-top: 1px solid lightgrey;"),
                      inputDataVisual
                  )
                  
)
## 2. page 2 #########################
## 2.1 parameter conrtrol panel ######
## ++++ 2.1.1 CV parameter control ########
parameterBox1 <- fluidRow(
  div(align = "middle",h3('Model Parameters for G2P and Cross Validation'),width = "50%"),
  column(12,radioButtons("GPorCV", h4("Genotype to phenotype prediction (G2P) or CrossValidation (CV):"), 
               choices = list("G2P" = 'G2P', "CV" = 'CV'),inline = T)),
  conditionalPanel(
    condition = "input.GPorCV == 'G2P'",
         column(4,
                h4('Input training sets index:'),
                textInput("trainIdx", "Training set index: ",  value = ""),
                h5('e.g.: 1:100, c(1:10,40:100)')),
         column(4,
                h4('Input testing sets index:'),
                textInput("testIdx", "Testing set index: ",  value = ""),
                h5('e.g.: 101:200, c(101:110,140:200)')),
         column(4,
                h4('Select modeling trait:'),
                uiOutput("traitSelectionPanel1")
                )
  ),
  conditionalPanel(
    condition = "input.GPorCV == 'CV'",
    column(12,uiOutput("traitSelectionPanel2")),
    column(12,radioButtons("CVMethods", h4("CV methods: "), 
                           choices = list("Holdout" = 'holdout', "K-fold" = 'Kfold',"Leave-one-out" = 'LOOCV'),inline = T)),
    conditionalPanel(
      condition = "input.CVMethods == 'holdout'",
      column(6,numericInput("holdoutFold",label = "Proportion of testing-sets",value = 0.2,min = 0,max = 1)),
      column(6,numericInput("holdoutRepeat",label = "Repeats of CV",value = 10,min = 1,max = 9999))
    ),
    conditionalPanel(
      condition = "input.CVMethods == 'Kfold'",
      column(6,numericInput("KfoldFold",label = "Fold of K-fold CV",value = 5,min = 2,max = 50)),
      column(6,numericInput("KfoldRepeat",label = "Repeats of CV",value = 10,min = 1,max = 9999))
    )
  )
)
## ++++ 2.1.2 Models parameter control ########
parameterBox2 <- fluidRow(
  column(12,div(align = "left",h4('Setting of prediction model parameters'),width = "50%")),
  column(9,
         selectInput("modelsSelection", label = NULL,
                     list(`Statistics models` = c("BayesA", "BayesB", "BayesC","RRBLUP", "LASSO"),
                          `Machine-learning models` = c("RFR","SVR")
                     ),multiple  = T,selected = "RRBLUP",
         )),
  column(2,
   actionButton('showModelModal','Arguments')
  ),
  column(12,uiOutput("showSelectedModel")),
  column(12,conditionalPanel(
    condition = "input.selectedModels == 'BayesA' || input.selectedModels == 'BayesB' || input.selectedModels == 'BayesC'",
    h4("Models parameters for Bayes family:"),
    column(3,
           numericInput('G2P.nIter',label = "nIter",value = 1500)
    ),
    column(3,
           numericInput('G2P.burnIn',label = "burnIn",value = 500)
    ),
    column(3,
           numericInput('G2P.thin',label = 'thin',value = 5)
    ),
    column(3,
           numericInput('G2P.S0',label = 'S0',value = NULL)
    ),
    column(3,
           numericInput('G2P.df0',label = "df0",value = 5)
    ),
    column(3,
           numericInput('G2P.R2',label = "R2",value = 0.5)
    )
  )),
  column(12,conditionalPanel(
    condition = "input.selectedModels == 'RRBLUP'",
    HTML('<b style="color:red;"> The RR-BLUP model have none of extra parameter.</b>'),
  )),
  column(12,conditionalPanel(
    condition = "input.selectedModels == 'LASSO'",
    h4("Models parameters for LASSO:"),
    column(3,
           numericInput('G2P.alpha',label = 'alpha',value = 1)
    )
  )),
  column(12,conditionalPanel(
    condition = "input.selectedModels == 'SVR'",
    h4('Models parameters for SVR:'),
    column(3,
           numericInput('G2P.gamma',label = "gamma",value = 1)
    ),
    column(3,
           numericInput('G2P.cost',label = 'cost',value = 2^(-9))
    ),
    column(3,
           selectInput('G2P.kernel',"kernel",c("linear","polynomial","sigmoid","radial"))
    )
  )),
  column(12,conditionalPanel(
    condition = "input.selectedModels == 'RFR'",
    h4('Models parameters for RFR:'),
    column(3,
           numericInput('G2P.ntree',label = "ntree",value = 500)
    ),
    column(3,
           numericInput('G2P.nodesize',label = "nodesize",value = 1)
    )
  ))
)

## ++++ 2.1.3 Eval parameter control ########
parameterBox3 <- fluidRow(
  column(12,div(align = "left",h4('Setting of evaluation metric parameters'),width = "50%")),
  column(9,
         selectInput("measuresSelection", NULL,
              list(`Global` = c("Pearson", "Kendall", "Spearman","MSE", "R2"),
                   `Threshold-based` = c("RE","Kappa", "NDCG", "meanNDCG", "F-score","Accuracy")
                   ),multiple  = T,selected = "Pearson"
  )
  ),
  column(2,
         actionButton('showEvalModal','Arguments')
  ),
  column(12,uiOutput("showSelectedMeasures")),
  column(6,
       sliderInput('eval.topAlpha', 'Threshold: top alpha(%) ', min = 1, max = 100,value = c(1,90)),
  ),
  column(5,
       radioButtons('eval.BestIndividuals',"BestIndividuals",c("top","middle","bottom"),selected = "top",inline = T)
  ),
  column(12,
       h5("Threshold for global metrics?"),
       checkboxInput("eval.globalAlpha", label = "GlobalAlpha", value = FALSE)
  ),
  column(12,conditionalPanel(
    condition = "input.selectedMeasures == 'F-score'",
    column(6,
           numericInput('eval.Beta',label = 'Bata',value = 1)
    )
  ))
)

## ++++ 2.1.4 run button control #######
consoleBox <- fluidRow(
  div(align = "middle",h3('G2P console'),width = "50%"),
  br(),
  fluidRow(column(12,
         br(),
         column(2,
                div(align = "left",actionButton('G2P.run','G2P',icon("send outline icon"),class="btn btn-primary")),
                br()
         ),
         column(2,
                div(align = "left",actionButton('G2PCV.run','G2P-CV',icon("send outline icon"),class="btn btn-primary")),
                br()
         ),
         # br(),
         column(3,
                div(align = "left",actionButton('eval.run.G2P','EvaluationG2P',icon("send outline icon"),class="btn btn-primary"))
         ),
         column(3,
                div(align = "left",actionButton('eval.run.CV','EvaluationCV',icon("send outline icon"),class="btn btn-primary"))
         )
  )) 
)

## ++++ 2.2.1 selected arguments visulization ########
selectedArguments <- fluidRow(
  div(align = "middle",h3('Modeling, evaluation and cross validation parameters (current)'),width = "50%"),

  column(12,htmlOutput("selectedParametersGPorCV")),
  column(12,htmlOutput("dataInformations")),
  column(4,
         h4("Cross validation arguments:"),
         htmlOutput("selectedParametersCV")),
  column(4,
         h4("Modeling arguments:"),
         htmlOutput("selectedParametersModel")),
  column(4,
         h4("Evaluation arguments:"),
         htmlOutput("selectedParametersEval"))
)

# ++++ 2.1.4 visulization show ########
## ++++ 2.2.2 result table show ########
predResEvalRes <- tabsetPanel(type = "tabs",
                              tabPanel("Results data table",
                                       tabsetPanel(
                                       tabPanel("Datatable for G2P prediciton results",
                                       HTML('<b>Table:</b> The results of G2P<br/>'),
                                       br(),
                                       DTOutput('G2PPredRes') %>% withSpinner(color="#0dc5c1"),
                                       ## define box options )
                                       ),
                                       tabPanel("Datatable for G2P evaluation results",
                                                HTML('<b>Table:</b> The results of evaluation<br/>'),
                                                tabsetPanel(
                                                  tabPanel("Overview of evaluation results",
                                                           uiOutput("showSelectAlpha"),
                                                           DTOutput('evalResOverall') %>% withSpinner(color="#0dc5c1")),
                                                  tabPanel("Pearson",
                                                           DTOutput('evalResPearson') %>% withSpinner(color="#0dc5c1")),
                                                  tabPanel("Kendall",
                                                           DTOutput('evalResKendall') %>% withSpinner(color="#0dc5c1")),
                                                  tabPanel("Spearman",
                                                           DTOutput('evalResSpearman') %>% withSpinner(color="#0dc5c1")),
                                                  tabPanel("MSE",
                                                           DTOutput('evalResMSE') %>% withSpinner(color="#0dc5c1")),
                                                  tabPanel("R2",
                                                           DTOutput('evalResR2') %>% withSpinner(color="#0dc5c1")),
                                                  tabPanel("RE",
                                                           DTOutput('evalResRE') %>% withSpinner(color="#0dc5c1")),
                                                  tabPanel("Kappa",
                                                           DTOutput('evalResKappa') %>% withSpinner(color="#0dc5c1")),
                                                  tabPanel("F-score",
                                                           DTOutput('evalResF1') %>% withSpinner(color="#0dc5c1")),
                                                  tabPanel("Accuracy",
                                                           DTOutput('evalResAccuracy') %>% withSpinner(color="#0dc5c1")),
                                                  tabPanel("NDCG",
                                                           DTOutput('evalResNDCG') %>% withSpinner(color="#0dc5c1")),
                                                  tabPanel("meanNDCG",
                                                           DTOutput('evalResMeanNDCG') %>% withSpinner(color="#0dc5c1"))
                                                )
                                       )),
                                       style = "height:580px; overflow-y: scroll;overflow-x: scroll;"),
                              tabPanel("Visulization of results",
                                       tabsetPanel(
                                                   tabPanel("Scatter plot",
                                                            plotlyOutput('visualization.scatter',height = "100%",width = "auto") %>% withSpinner(color="#0dc5c1")
                                                   )
                                       )
                                       )
)

## ++++ 2.2.3 visulization control ########
visulizationControl <-  fluidRow(
  ##input isoform expression data
  h3('Select plot data:'),
  # selectInput('plotData','Plot data',c("Prediction results","Global evaluation results","Threshold evaluation results")),
  # conditionalPanel(
    condition = "input.plotData == 'Prediction results'",
                     # wellPanel(
                       h4('Plot parameter setting:'),
                       fluidRow(
                         column(6,
                                # textInput('method1', label = 'Method 1: ', value = ' '),
                                uiOutput("scatter.selectM1"),
                                selectInput('scatter.showline',"Show the fitting line?",c("FALSE","TRUE")),
                                sliderInput("scatter.alpha", "Diaphaneity", min = 0, max = 1,value = 0.8,sep = 0.01)
                         ),
                         column(6,
                                # textInput('method2', label = 'Method 2: ', value = ' '),
                                uiOutput("scatter.selectM2"),
                                selectInput('scatterColor_size',"Change of gradient?",c("TRUE","FALSE")),
                                sliderInput("scatter.sizeRange", "Points size", sep = 0.05,min = 0, max = 20,value = c(2,5))
                         ),
                         conditionalPanel(
                           condition = "input.scatterColor_size == 'TRUE'",
                           column(12,
                                  sliderInput('scatter.col.low', label = 'Color for min value:', min = 0,max = 8000,value = 5000),
                                  sliderInput('scatter.col.mid', label = 'Color for median value:', min = 0,max = 8000,value = 0),
                                  sliderInput('scatter.col.high', label = 'Color for max value:', min = 0,max = 8000,value = 7000)
                           )
                         ),
                         conditionalPanel(
                           condition = "input.scatterColor_size == 'FALSE'",
                           column(12,
                                  sliderInput('scatter.color',label = 'Color for points:',min = 0,max = 8000,value = 5000)
                           )
                         )
                       ),
                     # )
    # Only show this panel if Custom is selected
    conditionalPanel(
      condition = "input.predPlotType == 'Desity'",

      selectInput("predres.select.method", "Select the methods", c("BayesA","BayesB","RFR","SVC"))
    )
)

page2 <- tabPanel("Genotype to phenotype [modeling,prediction and cross validation]",
                  sidebarLayout(
                    sidebarPanel(
                      tabsetPanel(type = "tabs",
                      tabPanel("Modeling console",
                      parameterBox1,
                      hr(style = "border-top: 1px solid lightgrey;"),
                      parameterBox2,
                      hr(style = "border-top: 1px solid lightgrey;"),
                      parameterBox3,
                      hr(style = "border-top: 1px solid lightgrey;"),
                      consoleBox),
                      tabPanel("Visulization console",
                               visulizationControl)

                      ),
                      width = 5
                    ),
                    mainPanel(
                      width = 7,
                      wellPanel(selectedArguments),
                      predResEvalRes
                    )
                  )
)

### Final navigation bar #################
UI <- shinyUI(navbarPage(HTML("Interactive platform for genomic selection <b>(IP4GS)</b>"),
                   ##Page 1
                   page1,
                   ##page 2
                   page2,inverse = T,collapsible = T
))
