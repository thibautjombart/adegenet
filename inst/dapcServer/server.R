library(shiny)
library(adegenet)


## DEFINE THE SERVER SIDE OF THE APPLICATION
shinyServer(function(input, output) {

    ## GET DYNAMIC ANNOTATION
    graphTitle <- reactive({
      paste(input$dataset, ": DAPC scatterplot, axes ", input$xax,"-", input$yax, sep="")
    })

    ## DEFINE CAPTION
    output$caption <- renderText({
        graphTitle()
    })


    ## GET DATA ##
    getData <- reactive({
        out <- NULL

        if(input$datatype=="expl"){
            if(input$dataset=="microbov") data("microbov", package="adegenet", envir=environment())
            if(input$dataset=="sim2pop") data("sim2pop", package="adegenet", envir=environment())
            if(input$dataset=="nancycats") data("nancycats", package="adegenet", envir=environment())
            out <- get(input$dataset)
        }

        if(input$datatype=="file" && !is.null(input$datafile)){
            ## need to rename input file
            oldName <- input$datafile$datapath
            extension <- .readExt(input$datafile$name)
            newName <- paste(input$datafile$datapath, extension, sep=".")
            file.rename(oldName, newName)

            ## treat different types of input
            if(extension %in% c("gtx","gen","dat","GTX","GEN","DAT")){
                out <- import2genind(newName)
            }

            if(extension %in% c("RData","Rdata","Rda","rda")){
                out <- get(load(newName))
            }

            if(extension %in% c("fasta","fa","fas","aln","FASTA","FA","FAS","ALN")){
                out <- DNAbin2genind(fasta2DNAbin(newName))
            }
        }
        return(out)
    })


    ## DYNAMIC UI COMPONENTS ##
    ## SELECTION OF PCA AXES
    output$npca <- renderUI({
        if(!is.null(x <- getData())) {
            nmax <- min(dim(x@tab))
            def <- min(10, nmax)
        } else {
            nmax <- 1000
            def <- 1
        }
        sliderInput("npca", "Number of PCA axes retained:", min=1, max=nmax, value=def,step=1)
    })

    ## SELECTION OF DA AXES
    output$nda <- renderUI({
        if(!is.null(x <- getData())) {
            nmax <- length(levels(pop(x)))-1
            def <- min(2, nmax)
        } else {
            nmax <- 100
            def <- 1
        }
        sliderInput("nda", "Number of DA axes retained:", min=1, max=nmax, value=def,step=1)
    })


    ## PERFORM THE DAPC ##
    getDapc <- reactive({
        out <- NULL
        ## eval(bquote(data(input$dataset))) # not sure why this isn't working
        x <- getData()
        npca <- nda <- 1
        if(!is.null(input$npca)) npca <- input$npca
        if(!is.null(input$nda)) nda <- input$nda
        if(!is.null(x)) out <- dapc(x, n.pca=npca, n.da=nda)
        return(out)
    })


    ## ## GET PLOT PARAM ##
    ## getPlotParam <- reactive({
    ##     col.pal <- get(input$col.pal)
    ##     return(list(col.pal=col.pal))
    ## })


    ## MAKE OUTPUT PLOT ##
    output$scatterplot <- renderPlot({
        dapc1 <- getDapc()
        if(!is.null(dapc1)){

            ## get colors
            K <- length(levels(dapc1$grp))
            myCol <- get(input$col.pal)(K)

            ## get screeplot info
            scree.pca <- ifelse(input$scree.pca=="none", FALSE, TRUE)
            scree.da <- ifelse(input$scree.da=="none", FALSE, TRUE)
            cellipse <- ifelse(input$ellipses, 1.5, 0)
            cstar <- ifelse(input$ellipses, 1, 0)
            scatter(dapc1, xax=input$xax, yax=input$yax, col=myCol,
                    scree.pca=scree.pca, scree.da=scree.da,
                    posi.pca=input$scree.pca, posi.da=input$scree.da,
                    cellipse=cellipse, cstar=cstar, mstree=input$mstree,
                    cex=input$pointsize, clabel=input$labelsize, solid=1-input$alpha)
        }
    })


    ## MAKE SUMMARY PLOT ##
    output$summary <- renderPrint({
        dapc1 <- getDapc()
        if(!is.null(dapc1)){
            summary(dapc1)
        }
    })

    ## MAKE COMPOPLOT ##
    output$compoplot <- renderPlot({
        dapc1 <- getDapc()
        if(!is.null(dapc1)){
            ## get colors
            K <- length(levels(dapc1$grp))
            myCol <- get(input$col.pal)(K)
            ##myCol <- transp(myCol, 1-input$alpha)
            compoplot(dapc1, col=myCol, lab=input$compo.lab, legend=input$compo.legend)

        }
    })

})
