#######
## dapc
########
dapc <- function (x, ...) UseMethod("dapc")

###################
## dapc.data.frame
###################
dapc.data.frame <- function(x, grp, n.pca=NULL, n.da=NULL,
                            center=TRUE, scale=FALSE, var.contrib=TRUE, pca.info=TRUE,
                            pca.select=c("nbEig","percVar"), perc.pca=NULL, ..., dudi=NULL){

    ## FIRST CHECKS
    grp <- as.factor(grp)
    if(length(grp) != nrow(x)) stop("Inconsistent length for grp")
    pca.select <- match.arg(pca.select)
    if(!is.null(perc.pca) & is.null(n.pca)) pca.select <- "percVar"
    if(is.null(perc.pca) & !is.null(n.pca)) pca.select <- "nbEig"
    if(!is.null(dudi) && !inherits(dudi, "dudi")) stop("dudi provided, but not of class 'dudi'")


    ## SOME GENERAL VARIABLES
    N <- nrow(x)
    REDUCEDIM <- is.null(dudi)

    if(REDUCEDIM){ # if no dudi provided
        ## PERFORM PCA ##
        maxRank <- min(dim(x))
        pcaX <- dudi.pca(x, center = center, scale = scale, scannf = FALSE, nf=maxRank)
    } else { # else use the provided dudi
        pcaX <- dudi
    }
    cumVar <- 100 * cumsum(pcaX$eig)/sum(pcaX$eig)

    if(!REDUCEDIM){
        myCol <- rep(c("black", "lightgrey"), c(ncol(pcaX$li),length(pcaX$eig)))
    } else {
        myCol <- "black"
    }

    ## select the number of retained PC for PCA
    if(is.null(n.pca) & pca.select=="nbEig"){
        plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA", col=myCol)
        cat("Choose the number PCs to retain (>=1): ")
        n.pca <- as.integer(readLines(n = 1))
    }

    if(is.null(perc.pca) & pca.select=="percVar"){
        plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA", col=myCol)
        cat("Choose the percentage of variance to retain (0-100): ")
        nperc.pca <- as.numeric(readLines(n = 1))
    }

    ## get n.pca from the % of variance to conserve
    if(!is.null(perc.pca)){
        n.pca <- min(which(cumVar >= perc.pca))
        if(perc.pca > 99.999) n.pca <- length(pcaX$eig)
        if(n.pca<1) n.pca <- 1
    }


    ## keep relevant PCs - stored in XU
    X.rank <- sum(pcaX$eig > 1e-14)
    n.pca <- min(X.rank, n.pca)
    if(n.pca >= N) n.pca <- N-1
    if(n.pca > N/3) warning("number of retained PCs of PCA may be too large (> N /3)\n results may be unstable ")
    n.pca <- round(n.pca)

    U <- pcaX$c1[, 1:n.pca, drop=FALSE] # principal axes
    rownames(U) <- colnames(x) # force to restore names
    XU <- pcaX$li[, 1:n.pca, drop=FALSE] # principal components
    XU.lambda <- sum(pcaX$eig[1:n.pca])/sum(pcaX$eig) # sum of retained eigenvalues
    names(U) <- paste("PCA-pa", 1:ncol(U), sep=".")
    names(XU) <- paste("PCA-pc", 1:ncol(XU), sep=".")


    ## PERFORM DA ##
    ldaX <- lda(XU, grp, tol=1e-30) # tol=1e-30 is a kludge, but a safe (?) one to avoid fancy rescaling by lda.default
    lda.dim <- sum(ldaX$svd^2 > 1e-10)
    ldaX$svd <- ldaX$svd[1:lda.dim]
    ldaX$scaling <- ldaX$scaling[,1:lda.dim,drop=FALSE]

    if(is.null(n.da)){
        barplot(ldaX$svd^2, xlab="Linear Discriminants", ylab="F-statistic", main="Discriminant analysis eigenvalues", col=heat.colors(length(levels(grp))) )
        cat("Choose the number discriminant functions to retain (>=1): ")
        n.da <- as.integer(readLines(n = 1))
    }

    ##n.da <- min(n.da, length(levels(grp))-1, n.pca) # can't be more than K-1 disc. func., or more than n.pca
    n.da <- round(min(n.da, lda.dim)) # can't be more than K-1 disc. func., or more than n.pca
    predX <- predict(ldaX, dimen=n.da)


    ## BUILD RESULT
    res <- list()
    res$n.pca <- n.pca
    res$n.da <- n.da
    res$tab <- XU
    res$grp <- grp
    res$var <- XU.lambda
    res$eig <- ldaX$svd^2
    res$loadings <- ldaX$scaling[, 1:n.da, drop=FALSE]
    res$means <- ldaX$means
    res$ind.coord <-predX$x
    res$grp.coord <- apply(res$ind.coord, 2, tapply, grp, mean)
    res$prior <- ldaX$prior
    res$posterior <- predX$posterior
    res$assign <- predX$class
    res$call <- match.call()


    ## optional: store loadings of variables
    if(pca.info){
        res$pca.loadings <- as.matrix(U)
        res$pca.cent <- pcaX$cent
        res$pca.norm <- pcaX$norm
        res$pca.eig <- pcaX$eig
    }

    ## optional: get loadings of variables
    if(var.contrib){
        res$var.contr <- as.matrix(U) %*% as.matrix(ldaX$scaling[,1:n.da,drop=FALSE])
        f1 <- function(x){
            temp <- sum(x*x)
            if(temp < 1e-12) return(rep(0, length(x)))
            return(x*x / temp)
        }
        res$var.contr <- apply(res$var.contr, 2, f1)
    }

    class(res) <- "dapc"
    return(res)
} # end dapc.data.frame





#############
## dapc.matrix
#############
dapc.matrix <- function(x, ...){
    return(dapc(as.data.frame(x), ...))
}




#############
## dapc.genind
#############
dapc.genind <- function(x, pop=NULL, n.pca=NULL, n.da=NULL,
                        scale=FALSE, truenames=TRUE, var.contrib=TRUE, pca.info=TRUE,
                        pca.select=c("nbEig","percVar"), perc.pca=NULL, ...){

    ## FIRST CHECKS
    if(!is.genind(x)) stop("x must be a genind object.")

    if(is.null(pop)) {
        pop.fac <- pop(x)
    } else {
        pop.fac <- pop
    }

    if(is.null(pop.fac)) stop("x does not include pre-defined populations, and `pop' is not provided")


    ## SOME GENERAL VARIABLES
    N <- nrow(x@tab)

    ## PERFORM PCA ##
    maxRank <- min(dim(x@tab))

    X <- scaleGen(x, center = TRUE, scale = scale,
                  missing = "mean", truenames = truenames)

    ## CALL DATA.FRAME METHOD ##
    res <- dapc(X, grp=pop.fac, n.pca=n.pca, n.da=n.da,
                center=FALSE, scale=FALSE, var.contrib=var.contrib,
                pca.select=pca.select, perc.pca=perc.pca)

    res$call <- match.call()

    ## restore centring/scaling
    res$pca.cent <- attr(X, "scaled:center")

    if(scale) {
        res$pca.norm <- attr(X, "scaled:scale")
    }

    return(res)
} # end dapc.genind






######################
## Function dapc.dudi
######################
dapc.dudi <- function(x, grp, ...){
    return(dapc.data.frame(x$li, grp, dudi=x, ...))
}





#################
## dapc.genlight
#################
dapc.genlight <- function(x, pop=NULL, n.pca=NULL, n.da=NULL,
                          scale=FALSE,  var.contrib=TRUE, pca.info=TRUE,
                          pca.select=c("nbEig","percVar"), perc.pca=NULL, glPca=NULL, ...){
    ## FIRST CHECKS ##
    if(!inherits(x, "genlight")) stop("x must be a genlight object.")

    pca.select <- match.arg(pca.select)

    if(is.null(pop)) {
        pop.fac <- pop(x)
    } else {
        pop.fac <- pop
    }

    if(is.null(pop.fac)) stop("x does not include pre-defined populations, and `pop' is not provided")



    ## PERFORM PCA ##
    REDUCEDIM <- is.null(glPca)

    if(REDUCEDIM){ # if no glPca provided
        maxRank <- min(c(nInd(x), nLoc(x)))
        pcaX <- glPca(x, center = TRUE, scale = scale, nf=maxRank, loadings=FALSE, returnDotProd = TRUE, ...)
    }

    if(!REDUCEDIM){ # else use the provided glPca object
        if(is.null(glPca$loadings) & var.contrib) {
            warning("Contribution of variables requested but glPca object provided without loadings.")
            var.contrib <- FALSE
        }
        pcaX <- glPca
    }

    if(is.null(n.pca)){
        cumVar <- 100 * cumsum(pcaX$eig)/sum(pcaX$eig)
    }


    ## select the number of retained PC for PCA
    if(!REDUCEDIM){
        myCol <- rep(c("black", "lightgrey"), c(ncol(pcaX$scores),length(pcaX$eig)))
    } else {
        myCol <- "black"
    }

    if(is.null(n.pca) & pca.select=="nbEig"){
        plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA", col=myCol)
        cat("Choose the number PCs to retain (>=1): ")
        n.pca <- as.integer(readLines(n = 1))
    }

    if(is.null(perc.pca) & pca.select=="percVar"){
        plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA", col=myCol)
        cat("Choose the percentage of variance to retain (0-100): ")
        nperc.pca <- as.numeric(readLines(n = 1))
    }

    ## get n.pca from the % of variance to conserve
    if(!is.null(perc.pca)){
        n.pca <- min(which(cumVar >= perc.pca))
        if(perc.pca > 99.999) n.pca <- length(pcaX$eig)
        if(n.pca<1) n.pca <- 1
    }

    if(!REDUCEDIM){
        if(n.pca > ncol(pcaX$scores)) {
            n.pca <- ncol(pcaX$scores)
        }
    }


    ## recompute PCA with loadings if needed
    if(REDUCEDIM){
        pcaX <- glPca(x, center = TRUE, scale = scale, nf=n.pca, loadings=var.contrib, matDotProd = pcaX$dotProd)
    }


    ## keep relevant PCs - stored in XU
    N <- nInd(x)
    X.rank <- sum(pcaX$eig > 1e-14)
    n.pca <- min(X.rank, n.pca)
    if(n.pca >= N) n.pca <- N-1
    if(n.pca > N/3) warning("number of retained PCs of PCA may be too large (> N /3)\n results may be unstable ")

    U <- pcaX$loadings[, 1:n.pca, drop=FALSE] # principal axes
    XU <- pcaX$scores[, 1:n.pca, drop=FALSE] # principal components
    XU.lambda <- sum(pcaX$eig[1:n.pca])/sum(pcaX$eig) # sum of retained eigenvalues
    names(U) <- paste("PCA-pa", 1:ncol(U), sep=".")
    names(XU) <- paste("PCA-pc", 1:ncol(XU), sep=".")


    ## PERFORM DA ##
    ldaX <- lda(XU, pop.fac, tol=1e-30) # tol=1e-30 is a kludge, but a safe (?) one to avoid fancy rescaling by lda.default
    lda.dim <- sum(ldaX$svd^2 > 1e-10)
    ldaX$svd <- ldaX$svd[1:lda.dim]
    ldaX$scaling <- ldaX$scaling[,1:lda.dim,drop=FALSE]

    if(is.null(n.da)){
        barplot(ldaX$svd^2, xlab="Linear Discriminants", ylab="F-statistic", main="Discriminant analysis eigenvalues", col=heat.colors(length(levels(pop.fac))) )
        cat("Choose the number discriminant functions to retain (>=1): ")
        n.da <- as.integer(readLines(n = 1))
    }

    n.da <- min(n.da, length(levels(pop.fac))-1, n.pca, sum(ldaX$svd>1e-10)) # can't be more than K-1 disc. func., or more than n.pca
    n.da <- round(n.da)
    predX <- predict(ldaX, dimen=n.da)


    ## BUILD RESULT
    res <- list()
    res$n.pca <- n.pca
    res$n.da <- n.da
    res$tab <- XU
    res$grp <- pop.fac
    res$var <- XU.lambda
    res$eig <- ldaX$svd^2
    res$loadings <- ldaX$scaling[, 1:n.da, drop=FALSE]
    res$means <- ldaX$means
    res$ind.coord <-predX$x
    res$grp.coord <- apply(res$ind.coord, 2, tapply, pop.fac, mean)
    res$prior <- ldaX$prior
    res$posterior <- predX$posterior
    res$assign <- predX$class
    res$call <- match.call()


    ## optional: store loadings of variables
    if(pca.info){
        res$pca.loadings <- as.matrix(U)
        res$pca.cent <- glMean(x,alleleAsUnit=FALSE)
        if(scale) {
            res$pca.norm <- sqrt(glVar(x,alleleAsUnit=FALSE))
        } else {
            res$pca.norm <- rep(1, nLoc(x))
        }
        res$pca.eig <- pcaX$eig
    }

    ## optional: get loadings of variables
    if(var.contrib){
        res$var.contr <- as.matrix(U) %*% as.matrix(ldaX$scaling[,1:n.da,drop=FALSE])
        f1 <- function(x){
            temp <- sum(x*x)
            if(temp < 1e-12) return(rep(0, length(x)))
            return(x*x / temp)
        }
        res$var.contr <- apply(res$var.contr, 2, f1)
    }

    class(res) <- "dapc"
    return(res)
} # end dapc.genlight






######################
# Function print.dapc
######################
print.dapc <- function(x, ...){
    cat("\t#################################################\n")
    cat("\t# Discriminant Analysis of Principal Components #\n")
    cat("\t#################################################\n")
    cat("class: ")
    cat(class(x))
    cat("\n$call: ")
    print(x$call)
    cat("\n$n.pca:", x$n.pca, "first PCs of PCA used")
    cat("\n$n.da:", x$n.da, "discriminant functions saved")
    cat("\n$var (proportion of conserved variance):", round(x$var,3))
    cat("\n\n$eig (eigenvalues): ")
    l0 <- sum(x$eig >= 0)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5)
        cat(" ...\n\n")

    ## vectors
    TABDIM <- 4
    if(!is.null(x$pca.loadings)){
        TABDIM <- TABDIM + 3
    }
    sumry <- array("", c(TABDIM, 3), list(1:TABDIM, c("vector", "length", "content")))
    sumry[1, ] <- c('$eig', length(x$eig),  'eigenvalues')
    sumry[2, ] <- c('$grp', length(x$grp), 'prior group assignment')
    sumry[3, ] <- c('$prior', length(x$prior), 'prior group probabilities')
    sumry[4, ] <- c('$assign', length(x$assign), 'posterior group assignment')
    if(!is.null(x$pca.loadings)){
        sumry[5, ] <- c('$pca.cent', length(x$pca.cent), 'centring vector of PCA')
        sumry[6, ] <- c('$pca.norm', length(x$pca.norm), 'scaling vector of PCA')
        sumry[7, ] <- c('$pca.eig', length(x$pca.eig), 'eigenvalues of PCA')
    }
    class(sumry) <- "table"
    print(sumry)

    ## data.frames
    cat("\n")
    TABDIM <- 6
    if(!is.null(x$pca.loadings)){
        TABDIM <- TABDIM + 1
    }
    if(!is.null(x$var.contr)){
        TABDIM <- TABDIM + 1
    }

    sumry <- array("", c(TABDIM, 4), list(1:TABDIM, c("data.frame", "nrow", "ncol", "content")))

    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "retained PCs of PCA")
    sumry[2, ] <- c("$means", nrow(x$means), ncol(x$means), "group means")
    sumry[3, ] <- c("$loadings", nrow(x$loadings), ncol(x$loadings), "loadings of variables")
    sumry[4, ] <- c("$ind.coord", nrow(x$ind.coord), ncol(x$ind.coord), "coordinates of individuals (principal components)")
    sumry[5, ] <- c("$grp.coord", nrow(x$grp.coord), ncol(x$grp.coord), "coordinates of groups")
    sumry[6, ] <- c("$posterior", nrow(x$posterior), ncol(x$posterior), "posterior membership probabilities")
    if(!is.null(x$pca.loadings)){
        sumry[7, ] <- c("$pca.loadings", nrow(x$pca.loadings), ncol(x$pca.loadings), "PCA loadings of original variables")
    }
    if(!is.null(x$var.contr)){
        sumry[TABDIM, ] <- c("$var.contr", nrow(x$var.contr), ncol(x$var.contr), "contribution of original variables")
    }
    class(sumry) <- "table"
    print(sumry)

    ## cat("\nother elements: ")
    ## if (length(names(x)) > 15)
    ##     cat(names(x)[15:(length(names(x)))], "\n")
    ## else cat("NULL\n")
    cat("\n")
} # end print.dapc






##############
## summary.dapc
##############
summary.dapc <- function(object, ...){
    x <- object
    res <- list()

    ## number of dimensions
    res$n.dim <- ncol(x$loadings)
    res$n.pop <- length(levels(x$grp))

    ## assignment success
    temp <- as.character(x$grp)==as.character(x$assign)
    res$assign.prop <- mean(temp)
    res$assign.per.pop <- tapply(temp, x$grp, mean)

    ## group sizes
    res$prior.grp.size <- table(x$grp)
    res$post.grp.size <- table(x$assign)

    return(res)
} # end summary.dapc






##############
## scatter.dapc
##############
scatter.dapc <- function(x, xax=1, yax=2, grp=x$grp, col=seasun(length(levels(grp))), pch=20, bg="white", solid=.7,
                         scree.da=TRUE, scree.pca=FALSE, posi.da="bottomright", posi.pca="bottomleft", bg.inset="white",
                         ratio.da=.25, ratio.pca=.25, inset.da=0.02, inset.pca=0.02, inset.solid=.5,
                         onedim.filled=TRUE, mstree=FALSE, lwd=1, lty=1, segcol="black",
                         legend=FALSE, posi.leg="topright", cleg=1, txt.leg=levels(grp),
                         cstar = 1, cellipse = 1.5, axesell = FALSE, label = levels(grp), clabel = 1, xlim = NULL, ylim = NULL,
                         grid = FALSE, addaxes = TRUE, origin = c(0,0), include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft",
                         cgrid = 1, pixmap = NULL, contour = NULL, area = NULL, ...){
    ONEDIM <- xax==yax | ncol(x$ind.coord)==1

    ## recycle color and pch
    col <- rep(col, length(levels(grp)))
    pch <- rep(pch, length(levels(grp)))
    col <- transp(col, solid)
    bg.inset <- transp(bg.inset, inset.solid)

    ## handle grp
    if(is.null(grp)){
        grp <- x$grp
    }

    if(!ONEDIM){
        ## set par
        opar <- par(mar = par("mar"))
        par(mar = c(0.1, 0.1, 0.1, 0.1), bg=bg)
        on.exit(par(opar))
        axes <- c(xax,yax)

        ## basic empty plot
        s.class(x$ind.coord[,axes], fac=grp, col=col, cpoint=0, cstar = cstar, cellipse = cellipse, axesell = axesell, label = label,
                clabel = clabel, xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, origin = origin, include.origin = include.origin,
                sub = sub, csub = csub, possub = possub, cgrid = cgrid, pixmap = pixmap, contour = contour, area = area)

        ## add points
        colfac <- pchfac <- grp
        levels(colfac) <- col
        levels(pchfac) <- pch
        colfac <- as.character(colfac)
        pchfac <- as.character(pchfac)
        if(is.numeric(col)) colfac <- as.numeric(colfac)
        if(is.numeric(pch)) pchfac <- as.numeric(pchfac)

        points(x$ind.coord[,xax], x$ind.coord[,yax], col=colfac, pch=pchfac, ...)
        s.class(x$ind.coord[,axes], fac=grp, col=col, cpoint=0, add.plot=TRUE, cstar = cstar, cellipse = cellipse, axesell = axesell, label = label,
                clabel = clabel, xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, origin = origin, include.origin = include.origin,
                sub = sub, csub = csub, possub = possub, cgrid = cgrid, pixmap = pixmap, contour = contour, area = area)

        ## add minimum spanning tree if needed
        if(mstree){
            meanposi <- apply(x$tab,2, tapply, grp, mean)
            D <- dist(meanposi)^2
            tre <- ade4::mstree(D)
            x0 <- x$grp.coord[tre[,1], axes[1]]
            y0 <- x$grp.coord[tre[,1], axes[2]]
            x1 <- x$grp.coord[tre[,2], axes[1]]
            y1 <- x$grp.coord[tre[,2], axes[2]]
            segments(x0, y0, x1, y1, lwd=lwd, lty=lty, col=segcol)
        }

    } else {
        ## set screeplot of DA to FALSE (just 1 bar)
        scree.da <- FALSE

        ## get plotted axis
        if(ncol(x$ind.coord)==1) {
            pcLab <- 1
        } else{
            pcLab <- xax
        }
        ## get densities
        ldens <- tapply(x$ind.coord[,pcLab], grp, density)
        allx <- unlist(lapply(ldens, function(e) e$x))
        ally <- unlist(lapply(ldens, function(e) e$y))
        par(bg=bg)
        plot(allx, ally, type="n", xlab=paste("Discriminant function", pcLab), ylab="Density")
        for(i in 1:length(ldens)){
            if(!onedim.filled) {
                lines(ldens[[i]]$x,ldens[[i]]$y, col=col[i], lwd=2) # add lines
            } else {
                polygon(c(ldens[[i]]$x,rev(ldens[[i]]$x)),c(ldens[[i]]$y,rep(0,length(ldens[[i]]$x))), col=col[i], lwd=2, border=col[i]) # add lines
            }
            points(x=x$ind.coord[grp==levels(grp)[i],pcLab], y=rep(0, sum(grp==levels(grp)[i])), pch="|", col=col[i]) # add points for indiv
        }
    }

    ## ADD INSETS ##
    ## group legend
    if(legend){
        ## add a legend
        temp <- list(...)$cex
        if(is.null(temp)) temp <- 1
        if(ONEDIM | temp<0.5 | all(pch=="")) {
            legend(posi.leg, fill=col, legend=txt.leg, cex=cleg, bg=bg.inset)
        } else {
            legend(posi.leg, col=col, legend=txt.leg, cex=cleg, bg=bg.inset, pch=pch, pt.cex=temp)
        }
    }

    ## eigenvalues discriminant analysis
    if(scree.da && ratio.da>.01) {
        inset <- function(){
            myCol <- rep("white", length(x$eig))
            myCol[1:x$n.da] <- "grey"
            myCol[c(xax, yax)] <- "black"
            myCol <- transp(myCol, inset.solid)
            barplot(x$eig, col=myCol, xaxt="n", yaxt="n", ylim=c(0, x$eig[1]*1.1))
            mtext(side=3, "DA eigenvalues", line=-1.2, adj=.8)
            box()
        }

        add.scatter(inset(), posi=posi.da, ratio=ratio.da, bg.col=bg.inset, inset=inset.da)
        ##add.scatter.eig(x$eig, ncol(x$loadings), axes[1], axes[2], posi=posi, ratio=ratio, csub=csub) # does not allow for bg
    }

    ## eigenvalues PCA
    if(scree.pca && !is.null(x$pca.eig) && ratio.pca>.01) {
        inset <- function(){
            temp <- 100* cumsum(x$pca.eig) / sum(x$pca.eig)
            myCol <- rep(c("black","grey"), c(x$n.pca, length(x$pca.eig)))
            myCol <- transp(myCol, inset.solid)
            plot(temp, col=myCol, ylim=c(0,115),
                 type="h", xaxt="n", yaxt="n", xlab="", ylab="", lwd=2)
            mtext(side=3, "PCA eigenvalues", line=-1.2, adj=.1)
        }
        add.scatter(inset(), posi=posi.pca, ratio=ratio.pca, bg.col=bg.inset, inset=inset.pca)
    }


    return(invisible(match.call()))
} # end scatter.dapc






############
## assignplot
############
assignplot <- function(x, only.grp=NULL, subset=NULL, new.pred=NULL, cex.lab=.75, pch=3){
    if(!inherits(x, "dapc")) stop("x is not a dapc object")

    ## handle data from predict.dapc ##
    if(!is.null(new.pred)){
        n.new <- length(new.pred$assign)
        x$grp <- c(as.character(x$grp), rep("unknown", n.new))
        x$assign <- c(as.character(x$assign), as.character(new.pred$assign))
        x$posterior <- rbind(x$posterior, new.pred$posterior)
    }


    ## treat other arguments ##
    if(!is.null(only.grp)){
        only.grp <- as.character(only.grp)
        ori.grp <- as.character(x$grp)
        x$grp <- x$grp[only.grp==ori.grp]
        x$assign <- x$assign[only.grp==ori.grp]
        x$posterior <- x$posterior[only.grp==ori.grp, , drop=FALSE]
    } else if(!is.null(subset)){
        x$grp <- x$grp[subset]
        x$assign <- x$assign[subset]
        x$posterior <- x$posterior[subset, , drop=FALSE]
    }


    ##table.paint(x$posterior, col.lab=ori.grp, ...)
    ## symbols(x$posterior)


    ## FIND PLOT PARAMETERS
    n.grp <- ncol(x$posterior)
    n.ind <- nrow(x$posterior)
    Z <- t(x$posterior)
    Z <- Z[,ncol(Z):1,drop=FALSE ]

    image(x=1:n.grp, y=seq(.5, by=1, le=n.ind), Z, col=rev(heat.colors(100)), yaxt="n", ylab="", xaxt="n", xlab="Clusters")
    axis(side=1, at=1:n.grp,tick=FALSE, labels=colnames(x$posterior))
    axis(side=2, at=seq(.5, by=1, le=n.ind), labels=rev(rownames(x$posterior)), las=1, cex.axis=cex.lab)
    abline(h=1:n.ind, col="lightgrey")
    abline(v=seq(0.5, by=1, le=n.grp))
    box()

    newGrp <- colnames(x$posterior)
    x.real.coord <- rev(match(x$grp, newGrp))
    y.real.coord <- seq(.5, by=1, le=n.ind)

    points(x.real.coord, y.real.coord, col="deepskyblue2", pch=pch)

    return(invisible(match.call()))
} # end assignplot





############
## compoplot
############
compoplot <- function(x, only.grp=NULL, subset=NULL, new.pred=NULL, col=NULL, lab=NULL,
                      legend=TRUE, txt.leg=NULL, ncol=4, posi=NULL, cleg=.8, bg=transp("white"), ...){
    if(!inherits(x, "dapc")) stop("x is not a dapc object")


    ## HANDLE ARGUMENTS ##
    ngrp <- length(levels(x$grp))

    ## col
    if(is.null(col)){
        col <- rainbow(ngrp)
    }

    ## lab
    if(is.null(lab)){
        lab <- rownames(x$tab)
    } else {
        ## recycle labels
       lab <- rep(lab, le=nrow(x$tab))
    }

    ## posi
    if(is.null(posi)){
        posi <- list(x=0, y=-.01)
    }

    ## txt.leg
    if(is.null(txt.leg)){
        txt.leg <- levels(x$grp)
    }

    ## HANDLE DATA FROM PREDICT.DAPC ##
    if(!is.null(new.pred)){
        n.new <- length(new.pred$assign)
        x$grp <- c(as.character(x$grp), rep("unknown", n.new))
        x$assign <- c(as.character(x$assign), as.character(new.pred$assign))
        x$posterior <- rbind(x$posterior, new.pred$posterior)
        lab <- c(lab, rownames(new.pred$posterior))
    }


    ## TREAT OTHER ARGUMENTS ##
    if(!is.null(only.grp)){
        only.grp <- as.character(only.grp)
        ori.grp <- as.character(x$grp)
        x$grp <- x$grp[only.grp==ori.grp]
        x$assign <- x$assign[only.grp==ori.grp]
        x$posterior <- x$posterior[only.grp==ori.grp, , drop=FALSE]
        lab <- lab[only.grp==ori.grp]
    } else if(!is.null(subset)){
        x$grp <- x$grp[subset]
        x$assign <- x$assign[subset]
        x$posterior <- x$posterior[subset, , drop=FALSE]
        lab <- lab[subset]
    }


    ## MAKE THE PLOT ##
    Z <- t(x$posterior)
    barplot(Z, border=NA, col=col, ylab="membership probability", names=lab, las=3, ...)

    if(legend){
        oxpd <- par("xpd")
        par(xpd=TRUE)
        legend(posi, fill=col, leg=txt.leg, cex=cleg, ncol=ncol, bg=bg)
        on.exit(par(xpd=oxpd))
    }

    return(invisible(match.call()))
} # end compoplot





###############
## a.score
###############
a.score <- function(x, n.sim=10, ...){
    if(!inherits(x,"dapc")) stop("x is not a dapc object")

    ## perform DAPC based on permuted groups
    lsim <- lapply(1:n.sim, function(i) summary(dapc(x$tab, sample(x$grp), n.pca=x$n.pca, n.da=x$n.da))$assign.per.pop)
    sumry <- summary(x)

    ## get the a-scores
    f1 <- function(Pt, Pf){
        tol <- 1e-7
        ##res <- (Pt-Pf) / (1-Pf)
        ##res[Pf > (1-tol)] <- 0
        res <- Pt-Pf
        return(res)
    }

    lscores <- lapply(lsim, function(e) f1(sumry$assign.per.pop, e))

    ## make a table of a-scores
    tab <- data.frame(lscores)
    colnames(tab) <- paste("sim", 1:n.sim, sep=".")
    rownames(tab) <- names(sumry$assign.per.pop)
    tab <- t(as.matrix(tab))

    ## make result
    res <- list()
    res$tab <- tab
    res$pop.score <- apply(tab, 2, mean)
    res$mean <- mean(tab)

    return(res)

} # end a.score







##############
## optim.a.score
##############
optim.a.score <- function(x, n.pca=1:ncol(x$tab), smart=TRUE, n=10, plot=TRUE,
                         n.sim=10, n.da=length(levels(x$grp)), ...){
    ## A FEW CHECKS ##
    if(!inherits(x,"dapc")) stop("x is not a dapc object")
    if(max(n.pca)>ncol(x$tab)) {
        n.pca <- min(n.pca):ncol(x$tab)
    }
    if(n.da>length(levels(x$grp))){
        n.da <- min(n.da):length(levels(x$grp))
    }
    pred <- NULL
    if(length(n.pca)==1){
        n.pca <- 1:n.pca
    }
    if(length(n.da)==1){
        n.da <- 1:n.da
    }


    ## AUXILIARY FUNCTION ##
    f1 <- function(ndim){
        temp <- dapc(x$tab[,1:ndim,drop=FALSE], x$grp, n.pca=ndim, n.da=x$n.da)
        a.score(temp, n.sim=n.sim)$pop.score
    }


    ## SMART: COMPUTE A FEW VALUES, PREDICT THE BEST PICK ##
    if(smart){
        if(!require(stats)) stop("the package stats is required for 'smart' option")
        o.min <- min(n.pca)
        o.max <- max(n.pca)
        n.pca <- pretty(n.pca, n) # get evenly spaced nb of retained PCs
        n.pca <- n.pca[n.pca>0 & n.pca<=ncol(x$tab)]
        if(!any(o.min==n.pca)) n.pca <- c(o.min, n.pca) # make sure range is OK
        if(!any(o.max==n.pca)) n.pca <- c(o.max, n.pca) # make sure range is OK
        lres <- lapply(n.pca, f1)
        names(lres) <- n.pca
        means <- sapply(lres, mean)
        sp1 <- smooth.spline(n.pca, means) # spline smoothing
        pred <- predict(sp1, x=1:max(n.pca))
        best <- pred$x[which.max(pred$y)]
    } else { ## DO NOT TRY TO BE SMART ##
        lres <- lapply(n.pca, f1)
        names(lres) <- n.pca
        best <- which.max(sapply(lres, mean))
        means <- sapply(lres, mean)
    }


    ## MAKE FINAL OUTPUT ##
    res <- list()
    res$pop.score <- lres
    res$mean <- means
    if(!is.null(pred)) res$pred <- pred
    res$best <- best

    ## PLOTTING (OPTIONAL) ##
    if(plot){
        if(smart){
            boxplot(lres, at=n.pca, col="gold", xlab="Number of retained PCs", ylab="a-score", xlim=range(n.pca)+c(-1,1), ylim=c(-.1,1.1))
            lines(pred, lwd=3)
            points(pred$x[best], pred$y[best], col="red", lwd=3)
            title("a-score optimisation - spline interpolation")
            mtext(paste("Optimal number of PCs:", res$best), side=3)
        } else {
            myCol <- rep("gold", length(lres))
            myCol[best] <- "red"
            boxplot(lres, at=n.pca, col=myCol, xlab="Number of retained PCs", ylab="a-score", xlim=range(n.pca)+c(-1,1), ylim=c(-.1,1.1))
            lines(n.pca, sapply(lres, mean), lwd=3, type="b")
            myCol <- rep("black", length(lres))
            myCol[best] <- "red"
            points(n.pca, res$mean, lwd=3, col=myCol)
            title("a-score optimisation - basic search")
            mtext(paste("Optimal number of PCs:", res$best), side=3)
        }
    }

    return(res)
} # end optim.a.score






#############
## as.lda.dapc
#############
as.lda <- function(...){
    UseMethod("as.lda")
}

as.lda.dapc <- function(x, ...){
    if(!inherits(x,"dapc")) stop("x is not a dapc object")
    res <- list()

    res$N <- nrow(res$ind.coord)
    res$call <- match.call()
    res$counts <- as.integer(table(x$grp))
    res$lev <- names(res$counts) <- levels(x$grp)
    res$means <- x$means
    res$prior <- x$prior
    res$scaling <- x$loadings
    res$svd <- sqrt(x$eig)

    class(res) <- "lda"

    return(res)
} # end as.lda.dapc






##############
## predict.dapc
##############
predict.dapc <- function(object, newdata, prior = object$prior, dimen,
                         method = c("plug-in", "predictive", "debiased"), ...){

    if(!inherits(object,"dapc")) stop("x is not a dapc object")
    method <- match.arg(method)

    x <- as.lda(object)


    ## HANDLE NEW DATA ##
    if(!missing(newdata)){
        ## make a few checks
        if(is.null(object$pca.loadings)) stop("DAPC object does not contain loadings of original variables. \nPlease re-run DAPC using 'pca.loadings=TRUE'.")
        newdata <- as.matrix(newdata) # to force conversion, notably from genlight objects
        if(ncol(newdata) != nrow(object$pca.loadings)) stop("Number of variables in newdata does not match original data.")

        ## centre/scale data
        for(i in 1:nrow(newdata)){ # this is faster for large, flat matrices)
            newdata[i,] <- (newdata[i,] - object$pca.cent) / object$pca.norm
        }
        newdata[is.na(newdata)] <- 0

        ## project as supplementary individuals
        XU <- newdata %*% as.matrix(object$pca.loadings)
    } else {
        XU <- object$tab
    }

    ## FORCE IDENTICAL VARIABLE NAMES ##
    colnames(XU) <- colnames(object$tab)


    ## HANDLE DIMEN ##
    if(!missing(dimen)){
        if(dimen > object$n.da) stop(paste("Too many dimensions requested. \nOnly", object$n.da, "discriminant functions were saved in DAPC."))
    } else {
        dimen <- object$n.da
    }

    ## CALL PREDICT.LDA ##
    temp <- predict(x, XU, prior, dimen, method, ...)


    ## FORMAT OUTPUT ##
    res <- list()
    res$assign <- temp$class
    res$posterior <- temp$posterior
    res$ind.scores <- temp$x

    return(res)

} # end predict.dapc







############
## xvalDapc
############

xvalDapc <- function (x, ...) UseMethod("xvalDapc")

xvalDapc.data.frame <- function(x, grp, n.pca.max, n.da=NULL, training.set = 0.9,
                                result=c("groupMean","overall"),
                                center=TRUE, scale=FALSE, n.pca=NULL, n.rep=10, ...){

    ## CHECKS ##
    grp <- factor(grp)
    n.pca <- n.pca[n.pca>0]
    result <- match.arg(result)
    if(is.null(n.da)) {
        n.da <- length(levels(grp))-1
    }

    ## GET TRAINING SET SIZE ##
    N <- nrow(x)
    N.training <- round(N*training.set)

    ## GET FULL PCA ##
    if(missing(n.pca.max)) n.pca.max <- min(dim(x))
    pcaX <- dudi.pca(x, nf=n.pca.max, scannf=FALSE, center=center, scale=scale)
    n.pca.max <- min(n.pca.max,pcaX$rank,N.training-1)

    ## DETERMINE N.PCA IF NEEDED ##
    if(is.null(n.pca)){
        n.pca <- round(pretty(1:n.pca.max,10))
    }
    n.pca <- n.pca[n.pca>0 & n.pca<(N.training-1)]

    ## FUNCTION GETTING THE % OF ACCURATE PREDICTION FOR ONE NUMBER OF PCA PCs ##
    ## n.pca is a number of retained PCA PCs
    VOID.GRP <- FALSE # will be TRUE if empty group happened
    get.prop.pred <- function(n.pca){
        f1 <- function(){
            toKeep <- sample(1:N, N.training)
            if(!(all(table(grp[toKeep])>0) & all(table(grp[-toKeep])>0))) VOID.GRP <<- TRUE
            temp.pca <- pcaX
            temp.pca$li <- temp.pca$li[toKeep,,drop=FALSE]
            temp.dapc <- suppressWarnings(dapc(x[toKeep,,drop=FALSE], grp[toKeep], n.pca=n.pca, n.da=n.da, dudi=temp.pca))
            temp.pred <- predict.dapc(temp.dapc, newdata=x[-toKeep,,drop=FALSE])
            if(result=="overall"){
                out <- mean(temp.pred$assign==grp[-toKeep])
            }
            if(result=="groupMean"){
                out <- mean(tapply(temp.pred$assign==grp[-toKeep], grp[-toKeep], mean), na.rm=TRUE)
            }
            return(out)
        }
        return(replicate(n.rep, f1()))
    }


    ## GET %SUCCESSFUL OF ACCURATE PREDICTION FOR ALL VALUES ##
    res.all <- unlist(lapply(n.pca, get.prop.pred))
    if(VOID.GRP) warning("At least one group was absent from the training / validating sets.\nTry using smaller training sets.")
    res <- data.frame(n.pca=rep(n.pca, each=n.rep), success=res.all)
    return(res)
} # end xvalDapc.data.frame


xvalDapc.matrix <- xvalDapc.data.frame






## #############
## ## discriVal
## #############

## discriVal <- function (x, ...) UseMethod("discriVal")

## discriVal.data.frame <- function(x, grp, n.pca.max, n.da=NULL, center=TRUE, scale=FALSE, n.pca=NULL, ...){

##     ## CHECKS ##
##     grp <- factor(grp)
##     n.pca <- n.pca[n.pca>0]
##     if(is.null(n.da)) {
##         n.da <- length(levels(grp))-1
##     }

##     ## GET FULL PCA ##
##     if(missing(n.pca.max)) n.pca.max <- min(dim(x))
##     pcaX <- dudi.pca(x, nf=n.pca.max, scannf=FALSE, center=center, scale=scale)
##     n.pca.max <- min(n.pca.max,pcaX$rank)

##     ## DETERMINE N.PCA IF NEEDED ##
##     if(is.null(n.pca)){
##         n.pca <- round(pretty(1:n.pca.max,10))
##     }
##     n.pca <- n.pca[n.pca>0 & n.pca<n.pca.max]

##     ## FUNCTION GETTING THE TOTAL DISCRIMINATION (SUM OF EIGENVALUES) FOR ONE GIVEN NB OF PCA PCs ##
##     ## n.pca is a number of retained PCA PCs
##     get.totdiscr <- function(n.pca){
##             temp.dapc <- suppressWarnings(dapc(x, grp, n.pca=n.pca, n.da=n.da, dudi=pcaX))
##             return(sum(temp.dapc$eig))
##     }


##     ## GET %SUCCESSFUL OF ACCURATE PREDICTION FOR ALL VALUES ##
##     res.all <- sapply(n.pca, get.totdiscr)
##     res <- data.frame(n.pca=n.pca, success=res.all)
##     return(res)
## } # end discriVal.data.frame


## discriVal.matrix <- discriVal.data.frame




## There's a bunch of problems down there, commenting it for nowé
## xval.dapc <- function(object, n.pca, n.da, training.set = 90, ...){
##   training.set = training.set/100
##   kept.id <- unlist(tapply(1:nInd(object), pop(object), function(e) {pop.size = length(e); pop.size.train = round(pop.size * training.set); sample(e, pop.size.train, replace=FALSE)})) # this can't work: nInd/pop not defined for DAPC objects
##   training <- object[kept.id]
##   validating <- object[-kept.id]
##   post = vector(mode = 'list', length = n.pca)
##   asgn = vector(mode = 'list', length = n.pca)
##   ind = vector(mode = 'list', length = n.pca)
##   mtch = vector(mode = 'list', length = n.pca)
##   for(i in 1:n.pca){
##     dapc.base = dapc(training, n.pca = i, n.da = 15) # Why 15??
##     dapc.p = predict.dapc(dapc.base, newdata = validating)
##     match.prp = mean(as.character(dapc.p$assign)==as.character(pop(validating)))
##     post[[i]] = dapc.p$posterior
##     asgn[[i]] = dapc.p$assign
##     ind[[i]] = dapc.p$ind.score
##     mtch[[i]] = match.prp
##   }
##   res = list(assign = asgn, posterior = post, ind.score = ind, match.prp = mtch)
##   return(res)
## } # end of xval.dapc

## xval.genind  <- function(object, n.pca, n.da, training.set = 90, ...){
##   res = xval.dapc(object = object, n.pca = n.pca, n.da = n.da, training.set = training.set)
##   return(res)
## }
## ###############
## ## randtest.dapc
## ###############
## ##randtest.dapc <- function(x, nperm = 999, ...){

## ##} # end randtest.dapc




######## TESTS IN R #######

## TEST PREDICT.DAPC ##
## data(sim2pop)
## temp <- seppop(sim2pop)
## temp <- lapply(temp, function(e) hybridize(e,e,n=30)) # force equal pop sizes
## hyb <- hybridize(temp[[1]], temp[[2]], n=30)
## newdat <- repool(temp[[1]], temp[[2]], hyb)
## pop(newdat) <- rep(c("pop A", "popB", "hyb AB"), c(30,30,30))


## ##dapc1 <- dapc(newdat[1:61],n.pca=10,n.da=1)
## dapc1 <- dapc(newdat[1:60],n.pca=2,n.da=1)
## scatter(dapc1)
## hyb.pred <- predict(dapc1, newdat[61:90])

## scatter(dapc1)
## points(hyb.pred$ind.scores, rep(.1, 30))

## assignplot(dapc1, new.pred=hyb.pred)
## title("30 indiv popA, 30 indiv pop B, 30 hybrids")
