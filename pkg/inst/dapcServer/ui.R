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
                                     conditionalPanel(
                                                      ## condition
                                                      condition = "$('li.active a').first().html()!='Help'",
                                                      radioButtons("datatype", "What data source to use?",
                                                                   list("Example from adegenet"="expl","Input file"="file"))),

                                     ## choice of dataset if source is an example
                                     conditionalPanel(condition = "input.datatype=='expl'&& $('li.active a').first().html()!= 'Help'",
                                                      selectInput("dataset", "Select an example dataset:",
                                                                  choices=c("microbov","sim2pop","nancycats"))
                                                      ),

                                     ## choice of dataset if source is a file
                                     conditionalPanel(condition = "input.datatype=='file'&& $('li.active a').first().html()!= 'Help'",
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
                                                      uiOutput("doxval"),
                                                      uiOutput("npcaMax")
                                                      ),


                                     ## Select Output variable:
                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()!='Cross-Validation'&& $('li.active a').first().html()!= 'Help'",
                                                      checkboxInput("useoptimnpca", "Use suggested number of PCA components?", FALSE)
                                                      ),


                                     ##sliderInput("npca", "Number of PCA axes retained:", min=1, max=1000, value=10),
                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()!='Cross-Validation'&& $('li.active a').first().html()!= 'Help'",
                                                      uiOutput("npca")
                                                      ),

                                     ## select number of DA axes
                                     ##sliderInput("nda", "Number of discriminant functions retained:", min=1, max=100, value=1),
                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()!= 'Help'",
                                                      uiOutput("nda")),

                                     ## nrep slider
                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()=='Cross-Validation'",
                                                      sliderInput("nrep",
                                                                  "Number of replicates:",
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
                                                                  value = 0.9,
                                                                  step = 0.01)
                                                      ),


                                     ## result type
                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()=='Cross-Validation'",
                                                      radioButtons("result", "Assess by:",
                                                                   list("Group" = "groupMean",
                                                                        "Overall" = "overall"))
                                                      ),



                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()=='Scatterplot'",
                                                      h3("Graphical parameters")
                                                      ),



                                     ## inputs specific of scatterplot tab
                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()=='Scatterplot'",

                                                      ## select first axis to plot
                                                      ##numericInput("xax", "Indicate the x axis", value=1, min=1),
                                                      uiOutput("xax"),

                                                      ## select second axis to plot
                                                      ##numericInput("yax", "Indicate the y axis", value=1, min=1),
                                                      uiOutput("yax")
                                                      ),



                                     conditionalPanel(
                                                      ## condition
                                                      condition = "$('li.active a').first().html()=='Scatterplot'",
                                                      h3("Aesthetics")
                                                      ),

                                     conditionalPanel(
                                                      ## condition
                                                      condition = "$('li.active a').first().html()=='Compoplot'",
                                                      h3("Aesthetics")
                                                      ),



                                     ## select color palette
                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()=='Scatterplot'",
                                                      selectInput("col.pal", "Indicate a color palette to be used",
                                                                  choices=c("funky","spectral","seasun","azur","wasp"))
                                                      ),

                                     ## select color palette
                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()=='Compoplot'",
                                                      selectInput("col.pal", "Indicate a color palette to be used",
                                                                  choices=c("funky","spectral","seasun","azur","wasp"))
                                                      ),



                                     ## select transparency
                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()=='Scatterplot'
                                        && $('li.active a').first().html()=='Compoplot'",
                                                      sliderInput("alpha", "Choose transparency", min=0, max=1, step=0.05, value=0.5)
                                                      ),



                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()==='Scatterplot'",
                                                      ## symbol size
                                                      sliderInput("pointsize", "Size of the points", value=1, min=0, max=10, step=0.2),

                                                      ## label size
                                                      sliderInput("labelsize", "Size of the labels", value=1, min=0, max=10, step=0.2),

                                                      ## add screeplot of PCA?
                                                      selectInput("screepca", "Position of the PCA screeplot:",
                                                                  choices=c("None" = "none",
                                                                  "Bottom right" = "bottomright",
                                                                  "Bottom left" = "bottomleft",
                                                                  "Top right" = "topright",
                                                                  "Top left" = "topleft")),

                                                      ## add screeplot of DA?
                                                      selectInput("screeda", "Position of the DA screeplot:",
                                                                  choices=c("None" = "none",
                                                                  "Bottom right" = "bottomright",
                                                                  "Bottom left" = "bottomleft",
                                                                  "Top right" = "topright",
                                                                  "Top left" = "topleft")),

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
                                                      ),

                                     ## inputs specific of loadingplot tab
                                     conditionalPanel(
                                                      ## condition
                                                      "$('li.active a').first().html()=='Loading Plot'",

                                                      h3("Graphical parameters"),

                                                      ## select axis to plot
                                                      uiOutput("LPax"),

                                                      checkboxInput("threshold", "Display threshold?", TRUE),

                                                      selectInput("thresholdMethod", "Method for selecting threshold:",
                                                                  choices=c("Third quantile" = "quantile",
                                                                  "Complete linkage clustering" = "complete",
                                                                  "Single linkage clustering" = "single",
                                                                  "Average linkage clustering" = "average",
                                                                  "Centroid clustering" = "centroid",
                                                                  "McQuitty's similarity analysis" = "mcquitty",
                                                                  "Median clustering" = "median",
                                                                  "Ward's minimum variance method" = "ward")),

                                                      checkboxInput("FS", "Select and describe features above threshold", FALSE))


                                     ), # end sidebarPanel

                        ## MAIN PANEL
                        mainPanel(
                                  tabsetPanel(

                                              tabPanel("Scatterplot",plotOutput("scatterplot")),

                                              tabPanel("Summary", verbatimTextOutput("summary")),

                                              tabPanel("Compoplot", plotOutput("compoplot")),

                                              tabPanel("Loading Plot", plotOutput("loadingplot"),
                                                       conditionalPanel(condition = "input.FS==1",
                                                                        h3("Number of selected vs. unselected alleles"),
                                                                        verbatimTextOutput("FS1"),
                                                                        h3("List of selected alleles"),
                                                                        verbatimTextOutput("FS2"),
                                                                        h3("Names of selected alleles"),
                                                                        verbatimTextOutput("FS3"),
                                                                        h3("Contributions of selected alleles to discriminant axis"),
                                                                        verbatimTextOutput("FS4"))),

                                              tabPanel("Cross-Validation", plotOutput("xvalPlot"),
                                                       h3("Mean success by number of PCs"),
                                                       verbatimTextOutput("xvalResults3"),
                                                       h3("Number of PCs with highest mean"),
                                                       verbatimTextOutput("xvalResults4"),
                                                       h3("RMSE by number of PCs"),
                                                       verbatimTextOutput("xvalResults5"),
                                                       h3("Number of PCs with lowest RMSE"),
                                                       verbatimTextOutput("xvalResults6"),
                                                       h3("Cross-validation results"),
                                                       verbatimTextOutput("xvalResults1"),
                                                       h3("Median and CI for random chance"),
                                                       verbatimTextOutput("xvalResults2")),





##################
                                              ## HELP SECTION ##
##################

                                              tabPanel("Help",

                                                       ## OVERVIEW ##

                                                       h3("Overview"),
                                                       p("Welcome to the DAPC Server Help Section. Under each section heading below you can find a brief description and useful information
                  about the content of each tab on the Server. If you want to know more about a specific item or are seeking the definition of a term,
                  you may be more interested in the Glossary on the adjacent tab."),

                                                       p("The DAPC Server aims to provide a user-friendly interactive application of some of the functions contained in the R package adegenet.
                  More information about adegenet can be found on the adegenet website:"),

                                                       a("http://adegenet.r-forge.r-project.org/", href="http://adegenet.r-forge.r-project.org/", target="_blank"),

                                                       p(br(),"On the DAPC Server, users can explore the Discriminant Analysis of Principal Components (DAPC) method.
                  DAPC is a multivariate method that uses genetic data to describe the differences between pre-defined biological populations.
                  DAPC is a dimension reduction approach that generates synthetic variables composed of weighted combinations of the original variables
                  in a dataset to optimally capture between-group variation. DAPC uses Principal Components Analysis (PCA) as a prior step
                  to Discriminant Analysis (DA) to identify up to (K - 1) linearly uncorrelated discriminant axes that can optimally discriminate between
                  K groups of individuals. Unlike DA alone, DAPC is able to perform this procedure when the number of variables (alleles) greatly
                  exceeds the number of individuals, and also in the presence of correlations between variables."),

                                                       p("DAPC is a supervised method, which means that the clusters of individuals to be analysed must be pre-defined by the user.
                  In cases where individuals are not classified into groups, procedures (like that of the find.clusters function in R) can be
                  used to identify clusters so that a DAPC analysis can be carried out. These procedures are not included in the present version of the
                  DAPC Server, however, which requires all data to be imported with inherently defined populations."),

                                                       p("For any dataset containing a set of genetic variables for any number of individuals and a population grouping factor,
                  DAPC can be used to explore between-population differentiation,
                  estimate the probabilities of individual assignment to all possible groups,
                  and to examine the contribution of individual alleles to population structuring."),

                                                       p("Acceptable file types for input datasets include 'gtx/gen/dat/GTX/GEN/DAT/RData/Rdata/Rda/rda' and 'GENETIX/genepop/Fstat/R data'."),


                                                       ## EMAIL ##
                                                       h3(br(),"Ask your questions on the adegenet forum"),
                                                       p("Use the adegenet forum to ask all non-confidential questions: ", a("send an email", href="mailto:adegenet-forum@lists.r-forge.r-project.org", target="_top")),
                                                       p("Make sure to describe your problem clearly and to provide, whenever possible, a reproducible example for any reported error. Please also give your session info (copy and paste the content of the serverInfo tab at the end of your email)"),
                                                       
                                                       p("Note that this mailing list is moderated, and if not subscribed your post may be differed by a day or two. To subscribe to the mailing list, go to", a("this page", href="http://lists.r-forge.r-project.org/cgi-bin/mailman/listinfo/adegenet-forum", target="_blank")),


                                                       ## CITATION ##

                                                       h3(br(),"Citing the DAPC Server"),
                                                       h5("Citation for adegenet:"),
                                                       p("Jombart T.(2008) adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics 24: 1403-1405. doi:10.1093/bioinformatics/btn129", a("[link to paper]", href="http://bioinformatics.oxfordjournals.org/cgi/reprint/btn129?ijkey=6sqx5BTXCdYtBZz&keytype=ref", target="_blank")),
                                                       p("Jombart T. and Ahmed I. (2011) adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics. doi: 10.1093/bioinformatics/btr521"),

                                                       h5("Citation for the DAPC:"),
                                                       p("Jombart T, Devillard S and Balloux F (2010) Discriminant analysis of principal components: a new method for the analysis of genetically structured populations. BMC Genetics 11:94. doi:10.1186/1471-2156-11-94", a("[link to paper]", href="http://www.biomedcentral.com/1471-2156/11/94", target="_blank")),


                                                       ## SCATTER PLOT ##

                                                       h3(br(),"Scatter Plot"),
                                                       p("The Scatter Plot page provides a visual assessment of between-population differentiation.
                  Generated by applying the R function scatterplot to a dapc object, the output generated will appear in one of two forms.
                  If only one DA is retained (always the case if there are only 2 groups),
                  or both the x-axis and y-axis of the scatterplot are set to the same value, the output will display the
                  densities of individuals on the given discriminant function.
                  If more than one DA is retained and selected, the output will display individuals as dots and
                  groups as inertia ellipses, and will represent the relative position of each along the two
                  selected axes."),

                                                       p("The number of axes retained in both the PCA and DA steps of DAPC will have an impact on the
                  analysis and affect the scatter plot. By default, the number of DA axes retained is set at
                  the maximum of (K - 1) axes, where K is the number of groups. The default value of the number of
                  PCA axes is more arbitrarily defined, however, the 'Use suggested number of PCA components?'
                  tickbox provides the user with the option to use cross-validation to identify and select an optimal number
                  of PCs, where one exists. For more on this, see the section on cross-validation."),

                                                       p("There are a wide variety of graphical parameters for the DAPC scatterplot that can be customised
                  by the user. Those parameters that lack intuitive definition are described further in the Glossary."),


                                                       ## SUMMARY ##

                                                       h3(br(),"Summary"),
                                                       p("This page provides a summary of the dapc object."),
                                                       p("$n.dim' indicates the number of retained DAPC axes, which is affected by both the number of PCA
                  axes and DA axes retained."),
                                                       p("'$n.pop' indicates the number of groups or populations, which is defined by the dataset."),
                                                       p("'$assign.prop' indicates the proportion of overall correct assignment"),
                                                       p("'$assign.per.pop' indicates the proportions of successful reassignment (based on the discriminant
                  functions) of individuals to their original clusters. Large values indicate clear-cut clusters, while low
                  values suggest admixed groups."),
                                                       p("'$prior.grp.size' indicates prior group sizes."),
                                                       p("'$post.grp.size' indicates posterior group sizes."),




                                                       ## COMPOPLOT ##

                                                       h3(br(),"Compoplot"),
                                                       p("This page displays a compoplot, which is a bar plot showing the probabilities of assignment of
                  individuals to the different clusters. Individuals are plotted along the x-axis and membership probabilities are
                  plotted along the y-axis.From the compoplot, one can draw inferences about potential
                  admixture, and about the way in which the selection of PCA axes affects the stability of membership probabilities."),


                                                       ## LOADING PLOT ##

                                                       h3("Loading Plot"),
                                                       p("The Loading Plot page allows the user to examine how the original variables contribute to the
                  discriminant functions created by DAPC. Variables are plotted along the x-axis, and the contribution
                  of those variables to the DAPC is plotted in the y-axis."),
                                                       p("The side panel on the Loading Plot page provides the option of selecting a threshold above which variables are identified.
                  This can be useful simply for clarifying the image; hence, by default, only variables above the third quartile threshold are labelled.
                  A drop-down menu contains a variety of clustering methods that can also be used to set this threshold.
                  If desired, the user can choose to 'Select and describe features above the threshold'"),


                                                       ## CROSS-VALIDATION ##

                                                       h3(br(),"Cross-validation"),
                                                       p("When the 'Perform cross-validation?' box is ticked, this optimisation procedure will be carried
                  out on the Cross-validation page."),
                                                       p("Cross-validation is an optimisation procedure that is used in the context of DAPC to identify the number of principal components
                  that gives rise to the model with the highest predictive capacity. In cross-validation for DAPC, the data is divided into a training
                  set and a validation set (by default, comprising 90% and 10% of the data, respectively). The analysis is run on the training set with
                  variable numbers of PCs retained, and the degree to which the analysis is able to accurately predict the group membership of
                  excluded individuals (those in the validation set) is used to select the optimal number of PCs to retain. This procedure is replicated
                  with different random sub-samples a number of times specified by a slider on the side panel. In  the interest of computational time,
                  only 3 replicates are performed by default, though more replicates are recommended to achieve greater optimisation.
                  Success is calculated either by group (the default) or measured as overall success."),

                                                       p("A scatterplot of the results is displayed, showing the number of PCs retained on the x-axis and
                  success on the y-axis. Individual replicates appear as dots, and the density of points is displayed
                  in blue."),
                                                       p("Ideally, the data should fall in an arc, indicating an optimal point at its maximum where
                  the number of PCs retained leads to better predictive success than numbers of PCs either
                  above or below."),
                                                       p("Below the plot, a variety of summary statistics are provided.
                  Ultimately, it is the number of PCs associated with the lowest RMSE (root mean squared error,
                  see Glossary) which is selected if 'Use suggested number of PCA components?' is ticked.",br(),br(),br())



                                                       ), # end Help section




                                              ## GLOSSARY ##

                                              tabPanel("Glossary",
                                                       h3("Compoplot"),
                                                       p("A compoplot is a bar plot showing the probabilities of assignment of individuals to the different clusters."),
                                                       h3("Cross-validation"),
                                                       p("Cross-validation is an optimisation procedure that is used in the context of DAPC to identify the number of principal components
                  that gives rise to the model with the highest predictive capacity. In cross-validation for DAPC, the data is divided into a training
                  set and a validation set (by default, comprising 90% and 10% of the data, respectively). The analysis is run on the training set with
                  variable numbers of PCs retained, and the degree to which the analysis is able to accurately predict the group membership of
                  excluded individuals (those in the validation set) is used to select the optimal number of PCs to retain."),
                                                       p("Note: Performing cross-validation can substantially improve the results of DAPC; however, the amount of computational time
                  required increases with the size of the dataset in question and the number of replicates carried out."),
                                                       h3("DA"),
                                                       p("Discriminant Analysis (DA) is a procedure for optimally describing the differences between groups of individuals.
                  DAPC uses DA subsequent to Principal Component Analysis to maximise discrimination between groups in conditions
                  where DA alone would be inappropriate."),
                                                       h3("DA axis"),
                                                       p("DAPC uses Discriminant Analysis (DA) to describe the differences between K groups of individuals along
                  a maximum of (K - 1) Discriminant Analysis axes (DA axes)."),
                                                       h3("DAPC"),
                                                       p("Discriminant Analysis of Principal Components (DAPC) is a multivariate method that uses genetic data to describe the
                  differences between pre-defined biological populations. DAPC uses Principal Component Analysis as a prior step to
                  Discriminant Analysis to identify weighted linear combinations of the original variables which give rise to optimal
                  between-group discrimination."),
                                                       p("For more, see: Jombart T, Devillard S and Balloux F (2010) Discriminant analysis of principal components: a new method for the
                  analysis of genetically structured populations. BMC Genetics11:94. doi:10.1186/1471-2156-11-94."),
                                                       h3("Inertia ellipse"),
                                                       p("The inertia ellipses displayed optionally on the DAPC scatter plot provide graphical summaries of a cloud of points.
                  They are meant to give shape to groups of individuals, and do not necessarily represent a 95% confidence
                  interval on the position of the centroid of a cluster."),
                                                       h3("Loading"),
                                                       p("Loadings provide a measure of the contribution of each original variable to the discrimination between groups along a given
                  discriminant axis. A loading plot is used to visualise these loadings so that one can, for example, assess the weight of each
                  variable and identify those variables whose contributions exceed a threshold of interest."),
                                                       h3("Minimum spanning tree"),
                                                       p("A minimum spanning tree is a graph theoretical phenomenon in which all of the vertices (in our case, cluster centroids) of
                  a graph are connected into the tree that contains the shortest set of possible paths between vertices."),
                                                       h3("PCA"),
                                                       p("Principal Component Analysis (PCA) is a multivariate statistical method that generates linearly
                  uncorrelated principal components (PCs) composed of weighted linear combinations of the original variables
                  to represent overall variation in the data in a reduced space."),
                                                       h3("PC, PCA axis"),
                                                       p("Principal components (PCs, or PCA axes) are sets of linearly uncorrelated synthetic variables composed of
                  weighted linear combinations of original variables. PCs are used to describe multivariate phenomenon in a smaller (or equal)
                  number of dimensions. The 'optimal' number of PCs depends on the data. In DAPC, retaining too few PCs will cause useful
                  information to be excluded from the analysis (and hence discriminative power to be lost), while retaining too many PCs
                  may lead to problems of overfitting and uninterpretability."),
                                                       h3("RMSE"),
                                                       p("The root mean squared error (RMSE) is a measure of the error of an estimator. On the DAPC Server, RMSE is used in
                  cross-validation to assess the ability of the model, with variable numbers of principal components retained, to achieve
                  perfect prediction of individuals into the correct group."),
                                                       p("RMSE = sqrt((1/n)*sum(i=1 to n)(Yhat.i - Yi)^2), where Yhat is a vector of n predictions,
                  and Y is the vector of the true values."),
                                                       h3("Screeplot"),
                                                       p("A graphical representation of the variance, or cumulative variance, contained in the set of principal components
                  or discriminant functions. Components shaded in black represent those that have been retained in the analysis."),
                                                       h3("Training set"),
                                                       p("The training set, in cross-validation, is the set of individuals retained in the analysis.
                  The complement of the training set is the 'validation set', which contains the remaining individuals excluded from the analysis
                  who are used to test the performance of the model which varies as a function of the number of PCs retained.",br(),br(),br())
                                                       ), # end Glossary


                                              ## SERVER INFO ##

                                              tabPanel("System info", verbatimTextOutput("systeminfo"))
                                              ) # end tabsetPanel
                                  ) # end mainPanel
                        ) # end pageWithSidebar
        ) # end shinyUI
