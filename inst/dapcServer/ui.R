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

                                     ## select number of PCA axes
                                     ##sliderInput("npca", "Number of PCA axes retained:", min=1, max=1000, value=10),
                                     uiOutput("npca"),

                                     ## select number of DA axes
                                     ##sliderInput("nda", "Number of discriminant functions retained:", min=1, max=100, value=1),
                                     uiOutput("nda"),

                                     ## select color palette
                                     selectInput("col.pal", "Indicate a color palette to be used",
                                                 choices=c("funky","spectral","seasun","azur","wasp")),

                                     ## select transparency
                                     sliderInput("alpha", "Choose transparency", min=0, max=1, step=0.05, value=0.5),

                                     ## inputs specific of scatterplot tab
                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()==='Scatterplot'",

                                                      ## select first axis to plot
                                                      numericInput("xax", "Indicate the x axis", value=1, min=1),

                                                      ## select second axis to plot
                                                      numericInput("yax", "Indicate the y axis", value=1, min=1),

                                                      ## symbol size
                                                      sliderInput("pointsize", "Size of the points", value=1, min=0, max=10, step=0.2),

                                                      ## label size
                                                      sliderInput("labelsize", "Size of the labels", value=1, min=0, max=10, step=0.2),

                                                      ## add screeplot of PCA?
                                                      selectInput("scree.pca", "Position of the PCA screeplot:",
                                                                  choices=c("none","bottomright","bottomleft","topright","topleft")),

                                                      ## add screeplot of DA?
                                                      selectInput("scree.da", "Position of the DA screeplot:",
                                                                  choices=c("none","bottomright","bottomleft","topright","topleft")),

                                                      ## plot ellipses?
                                                      checkboxInput("ellipses", "Show inertia ellipses?", value=TRUE),

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

                                              tabPanel("Compoplot", plotOutput("compoplot"))
                                              ) # end tabsetPanel
                                  ) # end mainPanel
                        ) # end pageWithSidebar
        )
