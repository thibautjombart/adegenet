library(shiny)
##library(adegenet)


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
            if(input$dataset=="sim2pop") data(sim2pop)
            if(input$dataset=="microbov") data(microbov)
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

            if(extension %in% c("RData","Rdata","Rda","RDA")){
                out <- get(load(newName))
            }

            if(extension %in% c("fasta","fa","fas","aln","FASTA","FA","FAS","ALN")){
                out <- DNAbin2genind(fasta2DNAbin(newName))
            }
        }
        return(out)
    })


    ## PERFORM THE DAPC ##
    getDapc <- reactive({
        out <- NULL
        ## eval(bquote(data(input$dataset))) # not sure why this isn't working
        x <- getData()
        if(!is.null(x)) out <- dapc(x, n.pca=input$n.pca, n.da=input$n.da)
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

            scatter(dapc1, xax=input$xax, yax=input$yax, col=myCol,
                    scree.pca=scree.pca, scree.da=scree.da,
                    posi.pca=input$scree.pca, posi.da=input$scree.da)
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

            compoplot(dapc1, col=myCol, lab=input$compo.lab, legend=input$compo.legend)

        }
    })

})
