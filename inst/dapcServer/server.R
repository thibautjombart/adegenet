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
    
    
    ###################################################################
    # # Cross-validation
    
    XVALDAPC <- reactive({
      if(input$Scale.n.pca){
      x <- getData()
      mat <- as.matrix(na.replace(x, method="mean"))
      grp <- pop(x)
      n.pca.max <- input$n.pca.max
      result <- input$result
      n.rep <- input$n.rep
      training.set <- input$training.set
      out <- xvalDapc(mat, grp, n.pca.max=n.pca.max,
                        result=result, n.rep=n.rep, n.da=n.da, training.set=training.set)
      }
      else{
        out <- NULL
      }
      return(out)
    })
    
    
    # xvalPlot
    output$xvalPlot <- renderPlot({
      XVAL1 <- XVALDAPC()
      if(!is.null(XVAL1)){
      x <- getData()
      mat <- as.matrix(na.replace(x, method="mean"))
      grp <- pop(x)
      XVAL2 <- XVAL1[[1]]
      successV <-as.vector(XVAL2$success)
      random <- replicate(300, mean(tapply(sample(grp)==grp, grp, mean)))
      q.GRP <- quantile(random, c(0.025,0.5,0.975))
      smoothScatter(XVAL2$n.pca, successV, nrpoints=Inf, pch=20, col=transp("black"),
                    ylim=c(0,1), xlab="Number of PCA axes retained",
                    ylab="Proportion of successful outcome prediction", 
                    main="DAPC Cross-Validation")
      print(abline(h=q.GRP, lty=c(2,1,2)))
      }
    })
    
    # xval Output  
    output$xvalResults3 <-renderPrint({ 
      XVAL1 <- XVALDAPC()
      if(!is.null(XVAL1)){
      print(XVAL1[[3]])
      }
      
    })
    output$xvalResults4 <-renderPrint({ 
      XVAL1 <- XVALDAPC()
      if(!is.null(XVAL1)){
      print(XVAL1[[4]])
      }
      
    })
    output$xvalResults5 <-renderPrint({ 
      XVAL1 <- XVALDAPC()
      if(!is.null(XVAL1)){
      print(XVAL1[[5]])
      }
      
    })
    output$xvalResults6 <-renderPrint({ 
      XVAL1 <- XVALDAPC()
      if(!is.null(XVAL1)){
      print(XVAL1[[6]])
      }
      
    })
    output$xvalResults1 <-renderPrint({ 
      XVAL1 <- XVALDAPC()
      if(!is.null(XVAL1)){
      print(XVAL1[[1]])
      }
      
    })
    output$xvalResults2 <-renderPrint({ 
      XVAL1 <- XVALDAPC()
      if(!is.null(XVAL1)){
      print(XVAL1[[2]])
      }
      
    })
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ###################################################################
    
    ## PERFORM THE DAPC ##
    getDapc <- reactive({
      out <- NULL
      ## eval(bquote(data(input$dataset))) # not sure why this isn't working
      x <- getData()
      npca <- nda <- 1
      
      
      # n.pca determined by xval or slider?
      if(input$Scale.n.pca==TRUE){
        XVAL1 <- XVALDAPC()
        n.pca <- as.integer(XVAL1[[6]])
      }
      else{
        if(!is.null(input$npca)) npca <- input$npca
      }
      
      # n.da determined by xval or slider?
        XVAL1 <- XVALDAPC()
        dapc1 <- XVAL1[[7]]
        n.da <- dapc1$n.da

      
      if(!is.null(x)) out <- dapc(x, n.pca=npca, n.da=nda, parallel=FALSE)
      return(out)
    })
    
    
 



    ## GET PLOT PARAM ##
    getPlotParam <- reactive({
        col.pal <- get(input$col.pal)
        return(list(col.pal=col.pal))
    })


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
