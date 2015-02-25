#############
## GENERIC ##
#############


#' Genetic transitive graphs
#' 
#' These functions are under development. Please email the author before using
#' them for published work.\cr
#' 
#' The function \code{gengraph} generates graphs based on genetic distances, so
#' that pairs of entities (individuals or populations) are connected if and
#' only if they are distant by less than a given threshold distance. Graph
#' algorithms and classes from the \code{\link[igraph]{igraph}} package are
#' used.\cr
#' 
#' \code{gengraph} is a generic function with methods for the following types
#' of objects:\cr - \code{matrix} (only numeric data)\cr - \code{dist} \cr -
#' \code{\linkS4class{genind}} objects (genetic markers, individuals)\cr -
#' \code{\linkS4class{genpop}} objects (genetic markers, populations)\cr -
#' \code{\link[ape]{DNAbin}} objects (DNA sequences)
#' 
#' 
#' @aliases gengraph gengraph.default gengraph.matrix gengraph.dist
#' gengraph.genind gengraph.genpop gengraph.DNAbin
#' @param x a \code{matrix}, \code{dist}, \code{\linkS4class{genind}},
#' \code{\linkS4class{genpop}}, or \code{DNAbin} object. For \code{matrix} and
#' \code{dist}, the object represents pairwise (by default, Hamming) distances
#' between considered individuals.
#' @param cutoff a \code{numeric} value indicating the cutoff point, i.e. the
#' distance at which two entities are no longer connected in the garph produced
#' by the method.
#' @param ngrp an \code{integer} indicating the number of groups to be looked
#' for. A message is issued if this exact number could not be found.
#' @param computeAll a \code{logical} stating whether to investigate solutions
#' for every (integer) cutoff point; defaults to FALSE.
#' @param plot a \code{logical} indicating whether plots should be drawn;
#' defaults to TRUE; this operation can take time for large, highly-connected
#' graphs.
#' @param show.graph a \code{logical} indicating whether the found graph should
#' be drawn, only used in the interactive mode; this operation can take time
#' for large, highly-connected graphs; defaults to FALSE.
#' @param col.pal a color palette used to define group colors.
#' @param method an \code{integer} ranging from 1 to 6 indicating the type of
#' method to be used to derive a matrix of pairwise distances between
#' populations; values from 1 to 5 are passed to the function
#' \code{dist.genpop}; 6 corresponds to pairwise Fst; other values are not
#' supported.
#' @param truenames a logical indicating whether original labels should be used
#' for plotting (TRUE), as opposed to indices of sequences (FALSE).
#' @param nbreaks an integer indicating the number of breaks used by the
#' heuristic when seeking an exact number of groups.
#' @param \dots further arguments to be used by other functions; currently not
#' used.
#' @return The class \code{gengraph} is a list with the following
#' components:\cr \item{graph}{a graph of class \code{\link[igraph]{igraph}}.}
#' \item{clust}{a list containing group information: \code{$membership}: an
#' integer giving group membership; \code{$csize}: the size of each cluster;
#' \code{$no}: the number of clusters} \item{cutoff}{the value used as a cutoff
#' point} \item{col}{the color used to plot each group.}
#' @author Original idea by Anne Cori and Christophe Fraser.  Implementation by
#' Thibaut Jombart \email{t.jombart@@imperial.ac.uk}.
#' @seealso The \code{\link[igraph]{igraph}} package.
#' @examples
#' 
#' \dontrun{
#' dat <- haploGen()
#' res <- gengraph(dat$seq, ngrp=1)
#' plot(res$graph)
#' }
#' 
#' @export gengraph
gengraph <-  function (x, ...) UseMethod("gengraph")





#############
## DEFAULT ##
#############
gengraph.default <- function(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE, plot=TRUE,
                             show.graph=TRUE, col.pal=funky, truenames=TRUE, nbreaks=10, ...){
    stop(paste("No method for objects of class",class(x)))
} # end gengraph.default





############
## MATRIX ##
############
##
## this is the basic method
##
gengraph.matrix <- function(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE, plot=TRUE, show.graph=TRUE, col.pal=funky,


#' Restore true labels of an object
#' 
#' The function \code{truenames} returns some elements of an object
#' (\linkS4class{genind} or \linkS4class{genpop}) using true names (as opposed
#' to generic labels) for individuals, markers, alleles, and population.\cr
#' 
#' 
#' @name truenames
#' @aliases truenames truenames-methods truenames,ANY-method
#' truenames,genind-method truenames,genpop-method
#' @docType methods
#' @param x a \linkS4class{genind} or a \linkS4class{genpop} object
#' @return If x\$pop is empty (NULL), a matrix similar to the x\$tab slot but
#' with true labels.
#' 
#' If x\$pop exists, a list with this matrix (\$tab) and a population vector
#' with true names (\$pop).\cr
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords manip
#' @examples
#' 
#' data(microbov)
#' microbov
#' 
#' microbov$tab[1:5,1:5]
#' truenames(microbov)$tab[1:5,1:5]
#' 
#' @export truenames
                            truenames=TRUE, nbreaks=10, ...){
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
        cutToTry <- pretty(x,nbreaks)
        cutToTry <- cutToTry[cutToTry>0 & cutToTry<max(x)]
        if(length(cutToTry)==0) cutToTry <- 1
        tempRes <- lapply(cutToTry, function(i) gengraph.matrix(x,cutoff=i))
        temp <- sapply(tempRes,function(e) e$clust$no)
        if(!min(abs(temp-ngrp))<1) warning(paste("The exact number of groups was not found. Tried increasing nbreaks"))
        cutoff <- cutToTry[which.min(abs(temp-ngrp))]
        ## if(!any(temp<ngrp)) {
        ##     cutoff <- 1
        ## } else {
        ##     cutoff <- cutToTry[max(which(temp>ngrp))]
        ## }

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
gengraph.dist <- function(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE, plot=TRUE,
                          show.graph=TRUE, col.pal=funky, truenames=TRUE, nbreaks=10, ...){
    ## CHECKS ##
    ## if(!require("igraph")) stop("igraph is required")

    ## USE MATRIX METHOD ##
    res <- gengraph(as.matrix(x), cutoff=cutoff, ngrp=ngrp, computeAll=computeAll, plot=plot, show.graph=show.graph, col.pal=col.pal,
                    truenames=truenames, nbreaks=nbreaks, ...)
    return(res)
} # end gengraph.dist







############
## GENIND ##
############
gengraph.genind <- function(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE, plot=TRUE,
                            show.graph=TRUE, col.pal=funky, truenames=TRUE, nbreaks=10, ...){
    ## CHECKS ##
    ## if(!require("igraph")) stop("igraph is required")

    ## COMPUTE DISTANCES ##
    x$tab[is.na(x$tab)] <- 0
    D <- (1-propShared(x))*nLoc(x)*ploidy(x)

    ## USE MATRIX METHOD ##
    res <- gengraph(D, cutoff=cutoff, ngrp=ngrp, computeAll=computeAll, plot=plot, show.graph=show.graph, col.pal=col.pal,
                    truenames=truenames, nbreaks=nbreaks, ...)
    if(truenames){
        V(res$graph)$label <- indNames(x)
    }

    return(res)
} # end gengraph.genind








############
## GENPOP ##
############
gengraph.genpop <- function(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE, plot=TRUE, show.graph=TRUE,
                            col.pal=funky, method=1, truenames=TRUE, nbreaks=10, ...){ ## CHECKS ##
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
                    truenames=truenames, nbreaks=nbreaks, ...)
    if(truenames){
        V(res$graph)$label <- x@pop.names
    }

    return(res)
} # end gengraph.genpop







############
## DNABIN ##
############
gengraph.DNAbin <- function(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE, plot=TRUE, show.graph=TRUE, col.pal=funky,
                            truenames=TRUE, nbreaks=10, ...){
    ## CHECKS ##
    ## if(!require("igraph")) stop("igraph is required")
    ## if(!require("ape")) stop("ape is required")

    ## COMPUTE DISTANCES ##
    D <- as.matrix(round(dist.dna(x,model="raw", pairwise.deletion = TRUE)*ncol(x)))

    ## USE MATRIX METHOD ##
    res <- gengraph(D, cutoff=cutoff, ngrp=ngrp, computeAll=computeAll, plot=plot, show.graph=show.graph, col.pal=col.pal,
                    truenames=truenames, nbreaks=nbreaks, ...)
    return(res)
} # end gengraph.DNAbin

