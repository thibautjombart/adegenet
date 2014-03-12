library(shiny)
##library(adegenet)

## DEFINE UI ##
shinyUI(
        pageWithSidebar(
                        ##  TITLE
                        headerPanel("DAPC web server"),

                        ## SIDE PANEL CONTENT
                        sidebarPanel(
                                     ## define the type of input
                                     radioButtons("datatype", "What data source to use?",
                                                  list("Example from adegenet"="expl","Input file"="file")),

                                     ## choice of dataset if source is an example
                                     conditionalPanel(condition = "input.datatype=='expl'",
                                                      selectInput("dataset", "Select an example dataset:",
                                                                  choices=c("microbov","sim2pop","nancycats"))
                                                      ),

                                     ## choice of dataset if source is a file
                                     conditionalPanel(condition = "input.datatype=='file'",
                                                      fileInput('datafile', 'Choose input file',
                                                                accept=c('gtx/gen/dat/GTX/GEN/DAT/RData/Rdata/Rda/rda', 'GENETIX/genepop/Fstat/R data')),
                                                      tags$hr()
                                                      ),



                                     ## CROSS-VALIDATION
                                     # n.pca.max slider
                                     conditionalPanel(
                                       ## condition
                                       "$('li.active a').first().html()=='Cross-Validation'",
                                       h3("Cross-validation"),
                                       sliderInput("n.pca.max",
                                                   "Max. number of PCs:",
                                                   min = 1,
                                                   max = 500,
                                                   value = 200)
                                     ),

                                     ## nrep slider
                                     conditionalPanel(
                                       ## condition
                                       "$('li.active a').first().html()=='Cross-Validation'",
                                       sliderInput("nrep",
                                                   "Number of repetitions:",
                                                   min = 1,
                                                   max = 100,
                                                   value = 3)
                                     ),

                                     ## trainingset slider
                                     conditionalPanel(
                                       ## condition
                                       "$('li.active a').first().html()=='Cross-Validation'",
                                       sliderInput("trainingset",
                                                   "Training set size:",
                                                   min = 0.1,
                                                   max = 0.95,
                                                   value = 0.9)
                                     ),


                                     ## result type
                                     conditionalPanel(
                                       ## condition
                                       "$('li.active a').first().html()=='Cross-Validation'",
                                       radioButtons("result", "Assess by:",
                                                    list("Group" = "groupMean",
                                                         "Overall" = "overall"))
                                     ),

                                     ## Select Output variable:
                                       checkboxInput("useoptimnpca", "Use suggested number of PCA components?", FALSE)
                                     ,


                                     ##sliderInput("npca", "Number of PCA axes retained:", min=1, max=1000, value=10),
                                     uiOutput("npca"),

                                     ## select number of DA axes
                                     ##sliderInput("nda", "Number of discriminant functions retained:", min=1, max=100, value=1),
                                     uiOutput("nda"),


                                     ## select color palette
                                     conditionalPanel(
                                       ## condition
                                       "$('li.active a').first().html()!='Cross-Validation'",
                                       selectInput("col.pal", "Indicate a color palette to be used",
                                                   choices=c("funky","spectral","seasun","azur","wasp"))
                                     ),


                                     ## select transparency
                                     conditionalPanel(
                                       ## condition
                                       "$('li.active a').first().html()!='Cross-Validation'",
                                       sliderInput("alpha", "Choose transparency", min=0, max=1, step=0.05, value=0.5)
                                     ),

                                     ## inputs specific of scatterplot tab
                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()==='Scatterplot'",

                                                      ## select first axis to plot
                                                      ##numericInput("xax", "Indicate the x axis", value=1, min=1),
                                                      uiOutput("xax"),

                                                      ## select second axis to plot
                                                      ##numericInput("yax", "Indicate the y axis", value=1, min=1),
                                                      uiOutput("yax"),

                                                      ## symbol size
                                                      sliderInput("pointsize", "Size of the points", value=1, min=0, max=10, step=0.2),

                                                      ## label size
                                                      sliderInput("labelsize", "Size of the labels", value=1, min=0, max=10, step=0.2),

                                                      ## add screeplot of PCA?
                                                      selectInput("screepca", "Position of the PCA screeplot:",
                                                                  choices=c("none","bottomright","bottomleft","topright","topleft")),

                                                      ## add screeplot of DA?
                                                      selectInput("screeda", "Position of the DA screeplot:",
                                                                  choices=c("none","bottomright","bottomleft","topright","topleft")),

                                                      ## plot ellipses?
                                                      checkboxInput("ellipses", "Show inertia ellipses?", value=TRUE),

                                                      ## plot stars?
                                                      checkboxInput("stars", "Link points to their centre?", value=TRUE),

                                                      ## plot minimum spanning tree?
                                                      checkboxInput("mstree", "Show minimum spanning tree?", value=FALSE)

                                                      ),

                                     ## input specific of compoplot tab
                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()==='Compoplot'",

                                                      ## add legend?
                                                      checkboxInput("compo.legend", "Add a legend?", TRUE),

                                                      ## add labels?
                                                      checkboxInput("compo.lab", "Display labels?", FALSE)
                                                      )
                                     ), # end sidebarPanel

                        ## MAIN PANEL
                        mainPanel(
                                  tabsetPanel(

                                    tabPanel("Scatterplot",plotOutput("scatterplot")),

                                    tabPanel("Summary", verbatimTextOutput("summary")),

                                    tabPanel("Compoplot", plotOutput("compoplot")),

                                    tabPanel("Cross-Validation", plotOutput("xvalPlot"),
                                             h3("Mean success by n.pca"),
                                             verbatimTextOutput("xvalResults3"),
                                             h3("n.pca with highest mean"),
                                             verbatimTextOutput("xvalResults4"),
                                             h3("RMSE by n.pca"),
                                             verbatimTextOutput("xvalResults5"),
                                             h3("n.pca with lowest RMSE"),
                                             verbatimTextOutput("xvalResults6"),
                                             h3("Cross-validation results"),
                                             verbatimTextOutput("xvalResults1"),
                                             h3("Median and CI for random chance"),
                                             verbatimTextOutput("xvalResults2"))
                                              ) # end tabsetPanel
                                  ) # end mainPanel
                        ) # end pageWithSidebar
        )
