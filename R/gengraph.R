#############
## GENERIC ##
#############
gengraph <-  function (x, ...) UseMethod("gengraph")





#############
## DEFAULT ##
#############
gengraph.default <- function(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE, plot=TRUE, show.graph=TRUE, col.pal=funky, truenames=TRUE, ...){
    stop(paste("No method for objects of class",class(x)))
} # end gengraph.default





############
## MATRIX ##
############
##
## this is the basic method
##
gengraph.matrix <- function(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE, plot=TRUE, show.graph=TRUE, col.pal=funky,
                            truenames=TRUE, ...){
    ## CHECKS ##
    ## if(!require("igraph")) stop("igraph is required")

    ## IF COMPUTEALL IS TRUE ##
    if(computeAll){
        cutoffvec <- 1:max(x)
        res <- lapply(cutoffvec, function(i) gengraph.matrix(x, cutoff=i, computeAll=FALSE))
        temp <- sapply(res, function(e) e$clust$no)
        if(plot){
            plot(cutoffvec, temp, xlab="Cut-off Hamming distance chosen", ylab="Number of groups")
        }
        return(res)
    }


    ## INTERACTIVE MODE IF BOTH CUTOFF AND NGRP MISSING ##
    if(is.null(cutoff) & is.null(ngrp)){
        chooseAgain <- TRUE
        while (chooseAgain) {
            if(plot){
                hist(x, nclass=50, col="deepskyblue1",xlab="Hamming distance",ylab="Frequency",main="Distribution of frequences")
            }
            cat("\nPlease choose a cutoff distance:  ")
            ans <- NA
            while(is.null(ans) || is.na(ans)) suppressWarnings(ans <- as.numeric(readLines(n = 1)))
            if(plot){
                abline(v=ans,col="red",lty=2, lwd=2)
            }
            res <- gengraph.matrix(x, cutoff=ans, truenames=truenames)
            if(truenames){
                V(res$graph)$label <- rownames(x)
            }

            cat(paste("\nNumber of clusters found:  ", res$clust$no, sep=""))
            if(plot && show.graph) plot(res$graph)
            ans <- ""
            while(!ans %in% c("y","n")){
                cat("\nAre you satisfied with this solution? (yes:y / no:n): ")
                ans <- tolower(readLines(n = 1))
            }
            if(ans=="y") chooseAgain <- FALSE
        }
        return(res)
    }



    ## MAIN CASE: IF CUT-OFF POINT IS GIVEN ##
    if(!is.null(cutoff)){
        x[x>=cutoff] <- 0
        g <- graph.adjacency(x, mode="undirected", weighted=TRUE, diag=FALSE)
        clust <- clusters(g)
        V(g)$color <- col.pal(clust$no)[clust$membership]
        col <- col.pal(clust$no)[1:clust$no]
        names(col) <- 1:clust$no

        ## assign labels to vertices
        if(truenames){
            V(g)$label <- rownames(x)
        } else {
            V(g)$label <- 1:nrow(x)
        }

        ## assign labels to edges
        if(length(E(g))>0) {
            E(g)$label <- E(g)$weight
        }

        ## make result
        res <- list(graph=g, clust=clusters(g), cutoff=cutoff, col=col)

    } else { ## IF CUT-OFF POINT NEEDS TO BE FOUND ##
        if(ngrp>=nrow(x)) stop("ngrp is greater than or equal to the number of individuals")

        ## FIRST HAVE A LOOK AT A RANGE OF VALUES ##
        cutToTry <- pretty(x,10)
        cutToTry <- cutToTry[cutToTry>1 & cutToTry<nrow(x)]
        if(length(cutToTry)==0) cutToTry <- 1
        tempRes <- lapply(cutToTry, function(i) gengraph.matrix(x,cutoff=i))
        temp <- sapply(tempRes,function(e) e$clust$no)
        if(!any(temp<ngrp)) {
            cutoff <- 1
        } else {
            cutoff <- cutToTry[max(which(temp>ngrp))]
        }

        ## FIND THE LOWEST CUTOFF GIVING NGRP ##
        res <- gengraph.matrix(x,cutoff=cutoff)
        while(res$clust$no>ngrp){
            cutoff <- cutoff+1
            res <- gengraph.matrix(x,cutoff=cutoff)
        }

        if(res$clust$no != ngrp) cat("\nNote: the exact number of clusters could not be found.\n")
    }


    ## RETURN ##
    ## if(truenames){
    ##     V(res$graph)$label <- rownames(x)
    ## } else {
    ##     V(res$graph)$label <- 1:nrow(x)
    ## }

    return(res)

} # end gengraph.matrix







############
## GENIND ##
############
gengraph.dist <- function(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE, plot=TRUE, show.graph=TRUE, col.pal=funky, truenames=TRUE, ...){
    ## CHECKS ##
    ## if(!require("igraph")) stop("igraph is required")

    ## USE MATRIX METHOD ##
    res <- gengraph(as.matrix(x), cutoff=cutoff, ngrp=ngrp, computeAll=computeAll, plot=plot, show.graph=show.graph, col.pal=col.pal,
                    truenames=truenames, ...)
    return(res)
} # end gengraph.dist







############
## GENIND ##
############
gengraph.genind <- function(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE, plot=TRUE, show.graph=TRUE, col.pal=funky,
                            truenames=TRUE, ...){
    ## CHECKS ##
    ## if(!require("igraph")) stop("igraph is required")

    ## COMPUTE DISTANCES ##
    x$tab[is.na(x$tab)] <- 0
    D <- (1-propShared(x))*nLoc(x)*ploidy(x)

    ## USE MATRIX METHOD ##
    res <- gengraph(D, cutoff=cutoff, ngrp=ngrp, computeAll=computeAll, plot=plot, show.graph=show.graph, col.pal=col.pal,
                    truenames=truenames, ...)
    if(truenames){
        V(res$graph)$label <- indNames(x)
    }

    return(res)
} # end gengraph.genind








############
## GENPOP ##
############
gengraph.genpop <- function(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE, plot=TRUE, show.graph=TRUE, col.pal=funky, method=1,
                            truenames=TRUE, ...){
    ## CHECKS ##
    ## if(!require("igraph")) stop("igraph is required")

    ## COMPUTE DISTANCES ##
    x$tab[is.na(x$tab)] <- 0
    if(method==6){
        D <- as.matrix(pairwise.fst(x))
    } else {
        D <- as.matrix(dist.genpop(x, method=method))
    }

    ## USE MATRIX METHOD ##
    res <- gengraph(D, cutoff=cutoff, ngrp=ngrp, computeAll=computeAll, plot=plot, show.graph=show.graph, col.pal=col.pal,
                    truenames=truenames, ...)
    if(truenames){
        V(res$graph)$label <- x@pop.names
    }

    return(res)
} # end gengraph.genpop







############
## DNABIN ##
############
gengraph.DNAbin <- function(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE, plot=TRUE, show.graph=TRUE, col.pal=funky,
                            truenames=TRUE, ...){
    ## CHECKS ##
    ## if(!require("igraph")) stop("igraph is required")
    ## if(!require("ape")) stop("ape is required")

    ## COMPUTE DISTANCES ##
    D <- as.matrix(round(dist.dna(x,model="raw", pairwise.deletion = TRUE)*ncol(x)))

    ## USE MATRIX METHOD ##
    res <- gengraph(D, cutoff=cutoff, ngrp=ngrp, computeAll=computeAll, plot=plot, show.graph=show.graph, col.pal=col.pal,
                    truenames=truenames, ...)
    return(res)
} # end gengraph.DNAbin

