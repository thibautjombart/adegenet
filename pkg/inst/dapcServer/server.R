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

          if(input$useoptimnpca){
            xval1 <- xvaldapc()
            npca <- as.integer(xval1[[6]])
            def <- npca}

        } else {
            nmax <- 1000
            def <- 1
        }
        sliderInput("npca", "Number of PCA axes retained:", min=1, max=nmax, value=def,step=1)
    })

    ## SELECTION OF DA AXES
    output$nda <- renderUI({
        if(!is.null(x <- getData())) {
            nmax <- max(length(levels(pop(x)))-1,2)
            def <- length(levels(pop(x)))-1
        } else {
            nmax <- 100
            def <- 1
        }
        sliderInput("nda", "Number of DA axes retained:", min=1, max=nmax, value=def,step=1)
    })

    ## SELECTION OF PLOTTED AXES
    output$xax <- renderUI({
        if(!is.null(x <- getData())) {
            nmax <- min(dim(x@tab))
        } else {
            nmax <- 1000
        }
        numericInput("xax", "Indicate the x axis", value=1, min=1, max=nmax)
    })

   output$yax <- renderUI({
       def <- 1
       nda <- 1
       if(!is.null(input$nda)) nda <- input$nda
       if(!is.null(x <- getData())) {
            nmax <- min(dim(x@tab))
            if(nda>1 && length(levels(pop(x)))>1) def <- 2
        } else {
            nmax <- 1000
        }

        numericInput("yax", "Indicate the y axis", value=def, min=1, max=nmax)
    })

    ## CROSS-VALIDATION
    ## DYNAMIC TICKBOX (TICKED IF OPTIMNPCA CHOSEN)
    output$doxval <- renderUI({
        checkboxInput("doxval", "Perform cross validation (computer intensive)?", input$doxval)
    })

   ## DYNAMIC SLIDER FOR MAX NPCA SELECTION
   output$npcaMax <- renderUI({
     if(!is.null(x <- getData())) {
       nmax <- min(dim(x@tab))
       def <- nmax
     } else {
       nmax <- 1000
       def <- 1
     }
     sliderInput("npcaMax", "Maximum number of PCs:", min=1, max=nmax, value=def,step=1)
   })


    ## CROSS-VALIDATION FUNCTION
    xvaldapc <- reactive({
        doxval <- FALSE
        if(!is.null(input$doxval)) doxval <- input$doxval
        if(input$useoptimnpca || doxval){
            x <- getData()
            mat <- as.matrix(na.replace(x, method="mean"))
            grp <- pop(x)
            n.pca.max <- input$n.pca.max
            result <- input$result
            n.rep <- input$nrep
            nda <- 1
            if(!is.null(input$nda)) nda <- input$nda
            training.set <- input$trainingset
            npcaMax <- 1
            if(!is.null(input$npcaMax)) npcaMax <- input$npcaMax
            out <- xvalDapc(mat, grp, n.pca.max=npcaMax,
                            result=result, n.rep=n.rep, n.da=nda, training.set=training.set,
                            xval.plot=FALSE)
        }
        else{
            out <- NULL
        }
        return(out)
    })


    ## XVALPLOT
    output$xvalPlot <- renderPlot({
        xval1 <- xvaldapc()
        if(!is.null(xval1)){
            x <- getData()
            mat <- as.matrix(na.replace(x, method="mean"))
            grp <- pop(x)
            xval2 <- xval1[[1]]
            successV <-as.vector(xval2$success)
            random <- replicate(300, mean(tapply(sample(grp)==grp, grp, mean)))
            q.GRP <- quantile(random, c(0.025,0.5,0.975))
            smoothScatter(xval2$n.pca, successV, nrpoints=Inf, pch=20, col=transp("black"),
                          ylim=c(0,1), xlab="Number of PCA axes retained",
                          ylab="Proportion of successful outcome prediction",
                          main="DAPC Cross-Validation")
            print(abline(h=q.GRP, lty=c(2,1,2)))
        }
    })

    ## XVAL OUTPUT
    output$xvalResults1 <-renderPrint({
        xval1 <- xvaldapc()
        if(!is.null(xval1)){
            print(xval1[[1]])
        }

    })
    output$xvalResults2 <-renderPrint({
        xval1 <- xvaldapc()
        if(!is.null(xval1)){
            print(xval1[[2]])
        }

    })
   output$xvalResults3 <-renderPrint({
        xval1 <- xvaldapc()
        if(!is.null(xval1)){
            print(xval1[[3]])
        }
    })
    output$xvalResults4 <-renderPrint({
        xval1 <- xvaldapc()
        if(!is.null(xval1)){
            print(xval1[[4]])
        }

    })
    output$xvalResults5 <-renderPrint({
        xval1 <- xvaldapc()
        if(!is.null(xval1)){
            print(xval1[[5]])
        }

    })
    output$xvalResults6 <-renderPrint({
        xval1 <- xvaldapc()
        if(!is.null(xval1)){
            print(xval1[[6]])
        }

    })


    ## PERFORM THE DAPC ##
    getDapc <- reactive({
        out <- NULL
        x <- getData()
        npca <- nda <- 1

        ## n.pca determined by xval or slider?
        if(input$useoptimnpca){
            xval1 <- xvaldapc()
            npca <- as.integer(xval1[[6]])
        } else {
            if(!is.null(input$npca)) npca <- input$npca
        }
        if(!is.null(input$nda)) nda <- input$nda

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
            scree.pca <- ifelse(input$screepca=="none", FALSE, TRUE)
            scree.da <- ifelse(input$screeda=="none", FALSE, TRUE)
            cellipse <- ifelse(input$ellipses, 1.5, 0)
            cstar <- ifelse(input$stars, 1, 0)
            scatter(dapc1, xax=input$xax, yax=input$yax, col=myCol,
                    scree.pca=scree.pca, scree.da=scree.da,
                    posi.pca=input$screepca, posi.da=input$screeda,
                    cellipse=cellipse, cstar=cstar, mstree=input$mstree,
                    cex=input$pointsize, clabel=input$labelsize, solid=1-input$alpha)
        } else {
            NULL
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


   ## DYNAMIC SELECTION OF DISCRIMINANT AXIS FOR LOADING PLOT
   output$LPax <- renderUI({
     def <- 1
     nda <- 1
     nmax <- 2
          if(!is.null(x <- getData())) {
            if(!is.null(input$nda)) nda <- input$nda
            nmax <- nda
            if(!is.null(input$LPax)) def <- input$LPax
          }

     numericInput("LPax", "Select discriminant axis", value=def, min=1, max=nmax)
   })



   # REACTIVE THRESHOLD/SNP SELECTION FUNCTION if using snpzip-like method
   selector <- reactive({

     dimension <- 1

     dapc1 <- getDapc()
     if(!is.null(dapc1)){
       if(!is.null(input$thresholdMethod)) method <- input$thresholdMethod
       if(!is.null(input$LPaxis)) dimension <- input$LPaxis
       x <- getData()
       mat <- as.matrix(na.replace(x, method="mean", quiet=TRUE))
     }

     if(method=="quantile"){
       x <- dapc1$var.contr[,dimension]
       thresh <- quantile(x,0.75)
       maximus <- which(x > thresh)
       n.snp.selected <- length(maximus)
       sel.snps <- mat[,maximus]
     }
     else{
     z <- dapc1$var.contr[,dimension]
     xTotal <- dapc1$var.contr[,dimension]
     toto <- which(xTotal%in%tail(sort(xTotal), 2000))
     z <- sapply(toto, function(e) xTotal[e])

     D <- dist(z)
     clust <- hclust(D,method)
     pop <- factor(cutree(clust,k=2,h=NULL))
     m <- which.max(tapply(z,pop,mean))
     maximus <- which(pop==m)
     maximus <- as.vector(unlist(sapply(maximus, function(e) toto[e])))
     popvect <- as.vector(unclass(pop))
     n.snp.selected <- sum(popvect==m)
     sel.snps <- mat[,maximus]
     }

     selection <- c((ncol(mat)-ncol(mat[,-maximus])), ncol(mat[,-maximus]))
     resultat <- list(selection, maximus, dimnames(sel.snps)[[2]], dapc1$var.contr[maximus, dimension])
     return(resultat)
   })



   ## MAKE LOADINGPLOT ##
   output$loadingplot <- renderPlot({
     dapc1 <- getDapc()
     LPaxis <- 1
     if(!is.null(dapc1)){
       ## get loadings for LP axis
       LPaxis <- 1
       if(!is.null(input$LPax)) LPaxis <- input$LPax
       if(input$threshold){
         # if threshold is by quantile
         if(input$thresholdMethod=="quantile"){
           x <- dapc1$var.contr[,LPaxis]
           def <- quantile(x,0.75)
         }else{
           # if threshold is by clustering
           select <- selector()
           thresh <- select[[2]]
           def <- abs(dapc1$var.contr[thresh][(which.min(dapc1$var.contr[thresh]))])-0.000001}
       }
       else{
         def <- NULL}
       loadingplot(dapc1$var.contr[,LPaxis], threshold=def)
     }
   })


   ## FEATURE SELECTION OUTPUT
   output$FS1 <-renderPrint({
     if(input$FS){
     fs1 <- selector()
     if(!is.null(fs1)){
       print(fs1[[1]])
     }
     }
   })

   output$FS2 <-renderPrint({
     if(input$FS){
     fs1 <- selector()
     if(!is.null(fs1)){
       print(fs1[[2]])
     }
     }
   })

   output$FS3 <-renderPrint({
     if(input$FS){
     fs1 <- selector()
     if(!is.null(fs1)){
       print(fs1[[3]])
     }
     }
   })

   output$FS4 <-renderPrint({
     if(input$FS){
     fs1 <- selector()
     if(!is.null(fs1)){
       print(fs1[[4]])
     }
     }
   })

    ## RENDER SYSTEM INFO ##
    output$systeminfo <- .render.server.info()

}) # end shinyServer
