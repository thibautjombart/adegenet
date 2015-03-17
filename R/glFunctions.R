
#' Auxiliary functions for genlight objects
#'
#' These functions provide facilities for usual computations using
#' \linkS4class{genlight} objects. When ploidy varies across individuals, the
#' outputs of these functions depend on whether the information units are
#' individuals, or alleles within individuals (see details).
#'
#' These functions are:
#'
#' - \code{glSum}: computes the sum of the number of second allele in each SNP.
#'
#' - \code{glNA}: computes the number of missing values in each SNP.
#'
#' - \code{glMean}: computes the mean number of second allele in each SNP.
#'
#' - \code{glVar}: computes the variance of the number of second allele in each
#' SNP.
#'
#' - \code{glDotProd}: computes dot products between (possibly centred/scaled)
#' vectors of individuals - uses compiled C code - used by glPca.
#'
#' === On the unit of information ===
#'
#' In the cases where individuals can have different ploidy, computation of
#' sums, means, etc. of allelic data depends on what we consider as a unit of
#' information.
#'
#' To estimate e.g. allele frequencies, unit of information can be considered
#' as the allele, so that a diploid genotype contains two samples, a triploid
#' individual, three samples, etc. In such a case, all computations are done
#' directly on the number of alleles. This corresponds to \code{alleleAsUnit =
#' TRUE}.
#'
#' However, when the focus is put on studying differences/similarities between
#' individuals, the unit of information is the individual, and all genotypes
#' possess the same information no matter what their ploidy is. In this case,
#' computations are made after standardizing individual genotypes to relative
#' allele frequencies. This corresponds to \code{alleleAsUnit = FALSE}.
#'
#' Note that when all individuals have the same ploidy, this distinction does
#' not hold any more.
#'
#' @aliases glSum glNA glMean glVar glDotProd
#' @param x a \linkS4class{genlight} object
#' @param alleleAsUnit a logical indicating whether alleles are considered as
#' units (i.e., a diploid genotype equals two samples, a triploid, three, etc.)
#' or whether individuals are considered as units of information.
#' @param center a logical indicating whether SNPs should be centred to mean
#' zero.
#' @param scale a logical indicating whether SNPs should be scaled to unit
#' variance.
#' @param useC a logical indicating whether compiled C code should be used
#' (TRUE) or not (FALSE, default).
#' @param parallel a logical indicating whether multiple cores -if available-
#' should be used for the computations (TRUE, default), or not (FALSE);
#' requires the package \code{parallel} to be installed (see details); this
#' option cannot be used alongside useCoption.
#' @param n.cores if \code{parallel} is TRUE, the number of cores to be used in
#' the computations; if NULL, then the maximum number of cores available on the
#' computer is used.
#' @return A numeric vector containing the requested information.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso - \code{\linkS4class{genlight}}: class of object for storing
#' massive binary SNP data.
#'
#' - \code{\link{dapc}}: Discriminant Analysis of Principal Components.
#'
#' - \code{\link{glPca}}: PCA for \linkS4class{genlight} objects.
#'
#' - \code{\link{glSim}}: a simple simulator for \linkS4class{genlight}
#' objects.
#'
#' - \code{\link{glPlot}}: plotting \linkS4class{genlight} objects.
#' @keywords multivariate
#' @examples
#'
#' \dontrun{
#' x <- new("genlight", list(c(0,0,1,1,0), c(1,1,1,0,0,1), c(2,1,1,1,1,NA)))
#' x
#' as.matrix(x)
#' ploidy(x)
#'
#' ## compute statistics - allele as unit ##
#' glNA(x)
#' glSum(x)
#' glMean(x)
#'
#' ## compute statistics - individual as unit ##
#' glNA(x, FALSE)
#' glSum(x, FALSE)
#' glMean(x, FALSE)
#'
#' ## explanation: data are taken as relative frequencies
#' temp <- as.matrix(x)/ploidy(x)
#' apply(temp,2, function(e) sum(is.na(e))) # NAs
#' apply(temp,2,sum, na.rm=TRUE) # sum
#' apply(temp,2,mean, na.rm=TRUE) # mean
#' }
#'

##########
## glSum
##########
## compute col sums
## removing NAs
##
glSum <- function(x, alleleAsUnit=TRUE, useC=FALSE){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    if(useC){
        ## use ploidy (sum absolute frequencies)
        if(alleleAsUnit){
            vecbyte <- unlist(lapply(x@gen, function(e) e$snp))
            nbVec <- sapply(x@gen, function(e) length(e$snp))
            nbNa <- sapply(NA.posi(x), length)
            naPosi <- unlist(NA.posi(x))
            res <- .C("GLsumInt", vecbyte, nbVec, length(x@gen[[1]]@snp[[1]]), nbNa, naPosi, nInd(x), nLoc(x), ploidy(x),
                      integer(nLoc(x)), PACKAGE="adegenet")[[9]]
        } else {
            ## sum relative frequencies
            vecbyte <- unlist(lapply(x@gen, function(e) e$snp))
            nbVec <- sapply(x@gen, function(e) length(e$snp))
            nbNa <- sapply(NA.posi(x), length)
            naPosi <- unlist(NA.posi(x))
            res <- .C("GLsumFreq", vecbyte, nbVec, length(x@gen[[1]]@snp[[1]]), nbNa, naPosi, nInd(x), nLoc(x), ploidy(x),
                      double(nLoc(x)), PACKAGE="adegenet")[[9]]
        }
    } else {
        ## use ploidy (sum absolute frequencies)
        if(alleleAsUnit){
            res <- integer(nLoc(x))
            for(e in x@gen){
                temp <- as.integer(e)
                temp[is.na(temp)] <- 0L
                res <- res + temp
            }
        } else {
            ## sum relative frequencies
            res <- numeric(nLoc(x))
            myPloidy <- ploidy(x)
            for(i in 1:nInd(x)){
                temp <- as.integer(x@gen[[i]]) / myPloidy[i]
                temp[is.na(temp)] <- 0
                res <- res + temp
            }
        }

    }
    names(res) <- locNames(x)
    return(res)

} # glSum





##########
## glNA
##########
## counts NB of NAs per column
##
## if alleleAsUnit, then effective is the number of alleles sampled (sum(ploidy(x)))
## otherwise, effective is simply the number of individuals (
glNA <- function(x, alleleAsUnit=TRUE){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## DEFAULT, VECTOR-WISE PROCEDURE ##
    res <- integer(nLoc(x))
    temp <- NA.posi(x)

    ## NAs in allele sampling
    if(alleleAsUnit){
        for(i in 1:length(temp)){
            if(length(temp[[i]])>0){
                res[temp[[i]]] <- res[temp[[i]]] + ploidy(x)[i]
            }
        }
    } else { ## NAs amongst individuals
        for(e in temp){
            if(length(e)>0){
                res[e] <- res[e] + 1
            }
        }
    }

    names(res) <- locNames(x)
    return(res)

} # glNA





##########
## glMean
##########
## computes SNPs means
## takes NAs into account
##
glMean <- function(x, alleleAsUnit=TRUE){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## DEFAULT, VECTOR-WISE PROCEDURE ##
    if(alleleAsUnit){ # use alleles
        N <- sum(ploidy(x)) - glNA(x, alleleAsUnit=TRUE)
        res <- glSum(x, alleleAsUnit=TRUE)/N
    } else { # use relative frequencies of individuals
        N <- nInd(x) - glNA(x, alleleAsUnit=FALSE)
        res <- glSum(x, alleleAsUnit=FALSE)/N
    }

    names(res) <- locNames(x)
    return(res)

} # glMean





########
## glVar
########
## computes SNPs variances
## takes NAs into account
##
glVar <- function(x, alleleAsUnit=TRUE){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## DEFAULT, VECTOR-WISE PROCEDURE ##
    res <- numeric(nLoc(x))
    myPloidy <- ploidy(x)

    if(alleleAsUnit){ # use alleles
        N <- sum(ploidy(x)) - glNA(x, alleleAsUnit=TRUE)
        xbar <- glMean(x, alleleAsUnit=TRUE)
        for(i in 1:nInd(x)){
            temp <- (as.integer(x@gen[[i]])/myPloidy[i] - xbar)^2
            temp[is.na(temp)] <- 0
            res <- res + temp*myPloidy[i]
        }
        res <- res/N
    } else { # use relative frequencies of individuals
        N <- nInd(x) - glNA(x, alleleAsUnit=FALSE)
        xbar <- glMean(x, alleleAsUnit=FALSE)

        for(i in 1:nInd(x)){
            temp <- (as.integer(x@gen[[i]])/myPloidy[i] - xbar)^2
            temp[is.na(temp)] <- 0L
            res <- res + temp
        }
        res <- res/N
    }

    names(res) <- locNames(x)
    return(res)

} # glVar






#############
## glDotProd
############
## computes all pairs of dot products
## between centred/scaled vectors
## of SNPs
glDotProd <- function(x, center=FALSE, scale=FALSE, alleleAsUnit=FALSE,
                      parallel=require("parallel"), n.cores=NULL){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## SOME CHECKS ##
    if(parallel && !require(parallel)) stop("parallel package requested but not installed")
    if(parallel && is.null(n.cores)){
        n.cores <- parallel::detectCores()
    }


    ## STORE USEFUL INFO ##
    N <- nInd(x)
    ind.names <- indNames(x)


    if(!parallel){ # DO NOT USE MULTIPLE CORES
        ## GET INPUTS TO C PROCEDURE ##
        if(center){
            mu <- glMean(x,alleleAsUnit=alleleAsUnit)
        } else {
            mu <- rep(0, nLoc(x))
        }

        if(scale){
            s <- sqrt(glVar(x,alleleAsUnit=alleleAsUnit))
            if(any(s<1e-10)) {
                warning("Null variances have been detected; corresponding alleles won't be standardized.")
            }
        } else {
            s <- rep(1, nLoc(x))
        }

        vecbyte <- unlist(lapply(x@gen, function(e) e$snp))
        nbVec <- sapply(x@gen, function(e) length(e$snp))
        nbNa <- sapply(NA.posi(x), length)
        naPosi <- unlist(NA.posi(x))
        lowerTriSize <- (nInd(x)*(nInd(x)-1))/2
        resSize <- lowerTriSize + nInd(x)

        ## CALL C FUNCTION ##
        temp <- .C("GLdotProd", vecbyte, nbVec, length(x@gen[[1]]@snp[[1]]), nbNa, naPosi, nInd(x), nLoc(x), ploidy(x),
                   as.double(mu), as.double(s), as.integer(!alleleAsUnit), double(resSize), PACKAGE="adegenet")[[12]]
    } else { # USE MULTIPLE CORES
        x <- seploc(x, n.block = n.cores) # one block per core (x is now a list of genlight)
        temp <- list()
        i <- 0
        for(block in x){
            i <- i+1
            ## GET INPUTS TO C PROCEDURE ##
            if(center){
                mu <- glMean(block,alleleAsUnit=alleleAsUnit)
            } else {
                mu <- rep(0, nLoc(block))
            }

            if(scale){
                s <- sqrt(glVar(block,alleleAsUnit=alleleAsUnit))
                if(any(s<1e-10)) {
                    warning("Null variances have been detected; corresponding alleles won't be standardized.")
                }
            } else {
                s <- rep(1, nLoc(block))
            }

            vecbyte <- unlist(lapply(block@gen, function(e) e$snp))
            nbVec <- sapply(block@gen, function(e) length(e$snp))
            nbNa <- sapply(NA.posi(block), length)
            naPosi <- unlist(NA.posi(block))
            lowerTriSize <- (nInd(block)*(nInd(block)-1))/2
            resSize <- lowerTriSize + nInd(block)

            ## CALL C FUNCTION ##
            temp[[i]] <- .C("GLdotProd", vecbyte, nbVec, length(block@gen[[1]]@snp[[1]]), nbNa, naPosi, nInd(block), nLoc(block), ploidy(block),
                       as.double(mu), as.double(s), as.integer(!alleleAsUnit), double(resSize), PACKAGE="adegenet")[[12]]
        }


        ## POOL BLOCK RESULTS TOGETHER ##
        temp <- Reduce("+", temp)
    }

    res <- temp[1:lowerTriSize]
    attr(res,"Size") <- N
    attr(res,"Diag") <- FALSE
    attr(res,"Upper") <- FALSE
    class(res) <- "dist"
    res <- as.matrix(res)
    diag(res) <- temp[(lowerTriSize+1):length(temp)]

    colnames(res) <- rownames(res) <- ind.names
    return(res)
} # end glDotProd






########
## glPca
########
##
## PCA for genlight objects
##


#' Principal Component Analysis for genlight objects
#' 
#' These functions implement Principal Component Analysis (PCA) for massive SNP
#' datasets stored as \linkS4class{genlight} object. This implementation has
#' the advantage of never representing to complete data matrix, therefore
#' making huge economies in terms of rapid access memory (RAM). When the
#' \code{parallel} package is available, \code{glPca} uses multiple-core
#' ressources for more efficient computations. \code{glPca} returns lists with
#' the class \code{glPca} (see 'value').
#' 
#' Other functions are defined for objects of this class:
#' 
#' - \code{print}: prints the content of a \code{glPca} object.
#' 
#' - \code{scatter}: produces scatterplots of principal components, with a
#' screeplot of eigenvalues as inset.
#' 
#' - \code{loadingplot}: plots the loadings of the analysis for one given axis,
#' using an adapted version of the generic function \code{loadingplot}.
#' 
#' === Using multiple cores ===
#' 
#' Most recent machines have one or several processors with multiple cores. R
#' processes usually use one single core. The package \code{parallel} allows
#' for parallelizing some computations on multiple cores, which can decrease
#' drastically computational time.
#' 
#' Lastly, note that using compiled C code (\code{useC=TRUE})is an alternative
#' for speeding up computations, but cannot be used together with the parallel
#' option.
#' 
#' @aliases glPca print.glPca scatter.glPca loadingplot.glPca
#' @param x for \code{glPca}, a \linkS4class{genlight} object; for
#' \code{print}, \code{scatter}, and \code{loadingplot}, a \code{glPca} object.
#' @param center a logical indicating whether the numbers of alleles should be
#' centered; defaults to TRUE
#' @param scale a logical indicating whether the numbers of alleles should be
#' scaled; defaults to FALSE
#' @param nf an integer indicating the number of principal components to be
#' retained; if NULL, a screeplot of eigenvalues will be displayed and the user
#' will be asked for a number of retained axes.
#' @param loadings a logical indicating whether loadings of the alleles should
#' be computed (TRUE, default), or not (FALSE). Vectors of loadings are not
#' always useful, and can take a large amount of RAM when millions of SNPs are
#' considered.
#' @param alleleAsUnit a logical indicating whether alleles are considered as
#' units (i.e., a diploid genotype equals two samples, a triploid, three, etc.)
#' or whether individuals are considered as units of information.
#' @param useC a logical indicating whether compiled C code should be used for
#' faster computations; this option cannot be used alongside parallel option.
#' @param parallel a logical indicating whether multiple cores -if available-
#' should be used for the computations (TRUE, default), or not (FALSE);
#' requires the package \code{parallel} to be installed (see details); this
#' option cannot be used alongside useCoption.
#' @param n.cores if \code{parallel} is TRUE, the number of cores to be used in
#' the computations; if NULL, then the maximum number of cores available on the
#' computer is used.
#' @param returnDotProd a logical indicating whether the matrix of dot products
#' between individuals should be returned (TRUE) or not (FALSE, default).
#' @param matDotProd an optional matrix of dot products between individuals,
#' NULL by default. This option is used internally to speed up computation time
#' when re-running the same PCA several times. Leave this argument as NULL
#' unless you really know what you are doing.
#' @param \dots further arguments to be passed to other functions.
#' @param xax,yax \code{integers} specifying which principal components should
#' be shown in x and y axes.
#' @param posi,bg,ratio arguments used to customize the inset in scatterplots
#' of \code{glPca} results. See \code{\link[ade4]{add.scatter}} documentation
#' in the ade4 package for more details.
#' @param
#' label,clabel,xlim,ylim,grid,addaxes,origin,include.origin,sub,csub,possub,cgrid,pixmap,contour,area
#' arguments passed to \code{\link[ade4]{s.class}}; see \code{?s.label} for
#' more information
#' @param at an optional numeric vector giving the abscissa at which loadings
#' are plotted. Useful when variates are SNPs with a known position in an
#' alignement.
#' @param threshold a threshold value above which values of x are identified.
#' By default, this is the third quartile of x.
#' @param axis an integer indicating the column of x to be plotted; used only
#' if x is a matrix-like object.
#' @param fac a factor defining groups of SNPs.
#' @param byfac a logical stating whether loadings should be averaged by groups
#' of SNPs, as defined by \code{fac}.
#' @param lab a character vector giving the labels used to annotate values
#' above the threshold.
#' @param cex.lab a numeric value indicating the size of annotations.
#' @param cex.fac a numeric value indicating the size of annotations for groups
#' of observations.
#' @param lab.jitter a numeric value indicating the factor of randomisation for
#' the position of annotations. Set to 0 (by default) implies no randomisation.
#' @param main the main title of the figure.
#' @param xlab the title of the x axis.
#' @param ylab the title of the y axis.
#' @param srt rotation of the labels; see ?text.
#' @param adj adjustment of the labels; see ?text.
#' @return === glPca objects ===
#' 
#' The class \code{glPca} is a list with the following components:\cr
#' \item{call}{the matched call.} \item{eig}{a numeric vector of eigenvalues.}
#' \item{scores}{a matrix of principal components, containing the coordinates
#' of each individual (in row) on each principal axis (in column).}
#' \item{loadings}{(optional) a matrix of loadings, containing the loadings of
#' each SNP (in row) for each principal axis (in column).}
#' 
#' -
#' 
#' === other outputs ===
#' 
#' Other functions have different outputs:\cr - \code{scatter} return the
#' matched call.\cr - \code{loadingplot} returns information about the most
#' contributing SNPs (see \code{\link{loadingplot.default}})
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso - \code{\linkS4class{genlight}}: class of object for storing
#' massive binary SNP data.
#' 
#' - \code{\link{glSim}}: a simple simulator for \linkS4class{genlight}
#' objects.
#' 
#' - \code{\link{glPlot}}: plotting \linkS4class{genlight} objects.
#' 
#' - \code{\link{dapc}}: Discriminant Analysis of Principal Components.
#' @keywords multivariate
#' @examples
#' 
#' \dontrun{
#' ## simulate a toy dataset
#' x <- glSim(50,4e3, 50, ploidy=2)
#' x
#' plot(x)
#' 
#' ## perform PCA
#' pca1 <- glPca(x, nf=2)
#' 
#' ## plot eigenvalues
#' barplot(pca1$eig, main="eigenvalues", col=heat.colors(length(pca1$eig)))
#' 
#' ## basic plot
#' scatter(pca1, ratio=.2)
#' 
#' ## plot showing groups
#' s.class(pca1$scores, pop(x), col=colors()[c(131,134)])
#' add.scatter.eig(pca1$eig,2,1,2)
#' }
#' 
#' @export glPca
glPca <- function(x, center=TRUE, scale=FALSE, nf=NULL, loadings=TRUE, alleleAsUnit=FALSE,
                  useC=TRUE, parallel=require("parallel"), n.cores=NULL,
                  returnDotProd=FALSE, matDotProd=NULL){
    if(!inherits(x, "genlight")) stop("x is not a genlight object")

    ## COMPUTE MEANS AND VARIANCES ##
    if(center) {
        vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
        if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
    }

    if(scale){
        vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
        if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
    }


    myPloidy <- ploidy(x)

    ## NEED TO COMPUTE DOT PRODUCTS ##
    if(is.null(matDotProd)){

        ## == if non-C code is used ==
        if(!useC){
            if(parallel && !require(parallel)) stop("parallel package requested but not installed")
            if(parallel && is.null(n.cores)){
                n.cores <- parallel::detectCores()
            }


            ## COMPUTE DOT PRODUCTS BETWEEN GENOTYPES ##
            ## to be fast, a particular function is defined for each case of centring/scaling

            ## NO CENTRING / NO SCALING
            if(!center & !scale){
                dotProd <- function(a,b, ploid.a, ploid.b){ # a and b are two SNPbin objects
                    a <- as.integer(a) / ploid.a
                    a[is.na(a)] <- 0
                    b <- as.integer(b) / ploid.b
                    b[is.na(b)] <- 0
                    return(sum( a*b, na.rm=TRUE))
                }
            }

            ## CENTRING / NO SCALING
            if(center & !scale){
                dotProd <- function(a,b, ploid.a, ploid.b){ # a and b are two SNPbin objects
                    a <- as.integer(a) / ploid.a
                    a[is.na(a)] <- vecMeans[is.na(a)]
                    b <- as.integer(b) / ploid.b
                    b[is.na(b)] <- vecMeans[is.na(b)]
                    return(sum( (a-vecMeans) * (b-vecMeans), na.rm=TRUE) )
                }
            }


            ## NO CENTRING / SCALING (odd option...)
            if(!center & scale){
                dotProd <- function(a,b, ploid.a, ploid.b){ # a and b are two SNPbin objects
                    a <- as.integer(a) / ploid.a
                    a[is.na(a)] <- 0
                    b <- as.integer(b) / ploid.b
                    b[is.na(b)] <- 0
                    return(sum( (a*b)/vecVar, na.rm=TRUE))
                }
            }


            ## CENTRING / SCALING
            if(center & scale){
                dotProd <- function(a,b, ploid.a, ploid.b){ # a and b are two SNPbin objects
                    a <- as.integer(a) / ploid.a
                    a[is.na(a)] <- vecMeans[is.na(a)]
                    b <- as.integer(b) / ploid.b
                    a[is.na(a)] <- vecMeans[is.na(a)]
                    return( sum( ((a-vecMeans)*(b-vecMeans))/vecVar, na.rm=TRUE ) )
                }
            }


            ## COMPUTE ALL POSSIBLE DOT PRODUCTS (XX^T / n) ##
            allComb <- combn(1:nInd(x), 2)
            if(parallel){
                allProd <- unlist(mclapply(1:ncol(allComb), function(i) dotProd(x@gen[[allComb[1,i]]], x@gen[[allComb[2,i]]], myPloidy[allComb[1,i]], myPloidy[allComb[2,i]]),
                                           mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE))
            } else {
                allProd <- unlist(lapply(1:ncol(allComb), function(i) dotProd(x@gen[[allComb[1,i]]], x@gen[[allComb[2,i]]], myPloidy[allComb[1,i]], myPloidy[allComb[2,i]]) ))
            }
            allProd <- allProd / nInd(x) # assume uniform weights

            ## shape result as a matrix
            attr(allProd,"Size") <- nInd(x)
            attr(allProd,"Diag") <- FALSE
            attr(allProd,"Upper") <- FALSE
            class(allProd) <- "dist"
            allProd <- as.matrix(allProd)

            ## compute the diagonal
            if(parallel){
                temp <- unlist(mclapply(1:nInd(x), function(i) dotProd(x@gen[[i]], x@gen[[i]], myPloidy[i], myPloidy[i]),
                                        mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE))/nInd(x)
            } else {
                temp <- unlist(lapply(1:nInd(x), function(i) dotProd(x@gen[[i]], x@gen[[i]], myPloidy[i], myPloidy[i]) ))/nInd(x)
            }
            diag(allProd) <- temp
        } else { # === use C computations ====
            allProd <- glDotProd(x, center=center, scale=scale, alleleAsUnit=alleleAsUnit, parallel=parallel, n.cores=n.cores)/nInd(x)
        }
    } else { # END NEED TO COMPUTE DOTPROD
        if(!all(dim(matDotProd)==nInd(x))) stop("matDotProd has wrong dimensions.")
        allProd <- matDotProd
    }

    ## PERFORM THE ANALYSIS ##
    ## eigenanalysis
    eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
    rank <- sum(eigRes$values > 1e-12)
    eigRes$values <- eigRes$values[1:rank]
    eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]

    ## scan nb of axes retained
    if(is.null(nf)){
        barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
        cat("Select the number of axes: ")
        nf <- as.integer(readLines(n = 1))
    }

    ## rescale PCs
    res <- list()
    res$eig <- eigRes$values
    nf <- min(nf, sum(res$eig>1e-10))
    ##res$matprod <- allProd # for debugging

    ## use: li = XQU = V\Lambda^(1/2)
    eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
    res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")


    ## GET LOADINGS ##
    ## need to decompose X^TDV into a sum of n matrices of dim p*r
    ## but only two such matrices are represented at a time
    if(loadings){
        if(scale) {
            vecSd <- sqrt(vecVar)
        }
        res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
        ## use: c1 = X^TDV
        ## and X^TV = A_1 + ... + A_n
        ## with A_k = X_[k-]^T v[k-]
        for(k in 1:nInd(x)){
            temp <- as.integer(x@gen[[k]]) / myPloidy[k]
            if(center) {
                temp[is.na(temp)] <- vecMeans[is.na(temp)]
                temp <- temp - vecMeans
            } else {
                temp[is.na(temp)] <- 0
            }
            if(scale){
                temp <- temp/vecSd
            }

            res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
        }

        res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
        res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
    }


    ## FORMAT OUTPUT ##
    colnames(res$scores) <- paste("PC", 1:nf, sep="")
    if(!is.null(indNames(x))){
        rownames(res$scores) <- indNames(x)
    } else {
         rownames(res$scores) <- 1:nInd(x)
    }

    if(!is.null(res$loadings)){
        colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
        if(!is.null(locNames(x)) & !is.null(alleles(x))){
            rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
        } else {
            rownames(res$loadings) <- 1:nLoc(x)
        }
    }

    if(returnDotProd){
        res$dotProd <- allProd
        rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
    }

    res$call <- match.call()

    class(res) <- "glPca"

    return(res)

} # glPca






###############
## print.glPca
###############
print.glPca <- function(x, ...){
    cat(" === PCA of genlight object ===")
    cat("\nClass: list of type glPca")
    cat("\nCall ($call):")
    print(x$call)
    cat("\nEigenvalues ($eig):\n", round(head(x$eig,6),3), ifelse(length(x$eig)>6, "...\n", "\n") )
    cat("\nPrincipal components ($scores):\n matrix with", nrow(x$scores), "rows (individuals) and", ncol(x$scores), "columns (axes)", "\n")
    if(!is.null(x$loadings)){
        cat("\nPrincipal axes ($loadings):\n matrix with", nrow(x$loadings), "rows (SNPs) and", ncol(x$loadings), "columns (axes)", "\n")
    }
    if(!is.null(x$dotProd)){
        cat("\nDot products between individuals ($dotProd):\n matrix with", nrow(x$dotProd), "rows and", ncol(x$dotProd), "columns", "\n")
    }
  cat("\n")
}





#################
## scatter.glPca
#################
scatter.glPca <- function(x, xax=1, yax=2, posi="bottomleft", bg="white", ratio=0.3,
                          label = rownames(x$scores), clabel = 1, xlim = NULL, ylim = NULL,
                          grid = TRUE, addaxes = TRUE, origin = c(0,0), include.origin = TRUE,
                          sub = "", csub = 1, possub = "bottomleft", cgrid = 1,
                          pixmap = NULL, contour = NULL, area = NULL, ...){

    ## set par
    opar <- par(mar = par("mar"))
    par(mar = c(0.1, 0.1, 0.1, 0.1), bg=bg)
    on.exit(par(opar))
    axes <- c(xax,yax)
    ## basic empty plot
    ## s.label(x$ind.coord[,axes], clab=0, cpoint=0, grid=FALSE, addaxes = FALSE, cgrid = 1, include.origin = FALSE, ...)
    s.label(x$scores[,axes], label = label, clabel = clabel, xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes,
            origin = origin, include.origin = include.origin, sub = sub, csub = csub, possub = possub, cgrid = cgrid,
            pixmap = pixmap, contour = contour, area = area)


    if(ratio>0.001) {
        add.scatter.eig(x$eig, ncol(x$scores), axes[1], axes[2], posi=posi, ratio=ratio, csub=csub)
    }

    return(invisible(match.call()))
} # end scatter.glPca






#####################
## loadingplot.glPca
#####################
loadingplot.glPca <- function(x, at=NULL, threshold=NULL, axis=1, fac=NULL, byfac=FALSE,
                        lab=rownames(x$loadings), cex.lab=0.7, cex.fac=1, lab.jitter=0,
                        main="Loading plot", xlab="SNP positions", ylab="Contributions", srt=90, adj=c(0,0.5), ... ){

    if(is.null(x$loadings)){
        warning("This object does not contain loadings. Re-run the analysis, specifying 'loadings=TRUE'.")
        return(invisible())
    }

    if(is.null(at)){
        at <- as.integer(gsub("[.]+.+$", "", rownames(x$loadings)))
    }
    if(is.null(threshold)){
        threshold <- quantile(x$loadings[,axis]^2,0.75)
    }

    res <- loadingplot.default(x$loadings^2, at=at, threshold=threshold, axis=axis, fac=fac, byfac=byfac,
                               lab=lab, cex.lab=cex.lab, cex.fac=cex.fac, lab.jitter=lab.jitter,
                               main=main, xlab=xlab, ylab=ylab, srt=srt, adj=adj, ...)

    axis(1)

    return(invisible(res))
} # end loadingplot.glPca





###############
## .glPca2dudi
###############
.glPca2dudi <- function(x){
    if(!inherits(x,"glPca")) stop("x is not a glPca object")
    old.names <- names(x)
    new.names <- sub("scores","li", old.names)
    new.names <- sub("loadings","c1", new.names)
    names(x) <- new.names
    class(x) <- c("pca","dudi")
    return(x)
} # end glPca2dudi




## TESTING ##
## x <- new("genlight", list(c(0,0,1,1,0), c(1,1,1,0,0,1), c(2,1,1,1,1,NA)))
## as.matrix(x)
## glNA(x)
## glSum(x)
## glMean(x)

## same ploidy everywhere
## x <- new("genlight", list(c(0,0,1,1,0), c(1,1,1,0,0,1), c(0,0,0,1,1,1)))
## f1 <- function(e) {return(mean((e-mean(e, na.rm=TRUE))^2, na.rm=TRUE))}
## all.equal(glVar(x), apply(as.matrix(x), 2, f1 )) # MUST BE TRUE
## all.equal(glVar(x,FALSE), apply(as.matrix(x), 2, f1 )) # MUST BE TRUE

## ## differences in ploidy
## x <- new("genlight", list(c(0,0,1,1,0), c(1,1,1,0,0,1), c(2,1,1,1,1,NA)))
## temp <- sweep(as.matrix(x), 1, c(1,1,2), "/")
## f2 <- function(e,w) {
##     mu <- weighted.mean(e, w, na.rm=TRUE)
##     res <- weighted.mean((e-mu)^2, w, na.rm=TRUE)
##     return(res)
## }

## all.equal(glMean(x), apply(temp,2,weighted.mean, w=c(1,1,2), na.rm=TRUE)) # MUST BE TRUE
## all.equal(glVar(x), apply(temp, 2, f2,w=c(1,1,2) )) # MUST BE TRUE

## all.equal(glMean(x,FALSE), apply(temp,2,mean,na.rm=TRUE)) # MUST BE TRUE
## all.equal(glVar(x,FALSE), apply(temp,2,f1)) # MUST BE TRUE



#### TESTING DOT PRODUCTS ####
## M <- matrix(sample(c(0,1), 100*1e3, replace=TRUE), nrow=100)
## rownames(M) <- paste("ind", 1:100)
## x <- new("genlight",M)

## ## not centred, not scaled
## res1 <- glDotProd(x,alleleAsUnit=FALSE, center=FALSE, scale=FALSE)
## res2 <- M %*% t(M)
## all.equal(res1,res2) # must be TRUE

## ##  centred, not scaled
## res1 <- glDotProd(x,alleleAsUnit=FALSE, center=TRUE, scale=FALSE)
## M <- scalewt(M,center=TRUE,scale=FALSE)
## res2 <- M %*% t(M)
## all.equal(res1,res2) # must be TRUE

## ##  centred, scaled
## res1 <- glDotProd(x,alleleAsUnit=FALSE, center=TRUE, scale=TRUE)
## M <- scalewt(M,center=TRUE,scale=TRUE)
## res2 <- M %*% t(M)
## all.equal(res1,res2) # must be TRUE


#### TESTING PCA ####
## M <- matrix(sample(c(0,1), 20*1000, replace=TRUE), nrow=20)
## rownames(M) <- paste("ind", 1:20)

## x <- new("genlight",M)
## res1 <- glPca(x, nf=4)
## res2 <- glPca(x, useC=FALSE, nf=4)
## res3 <- dudi.pca(M, center=TRUE,scale=FALSE, scannf=FALSE,nf=4)

## ## all must be TRUE
## all.equal(res1$eig,res3$eig)
## all.equal(res2$eig,res3$eig)
## all.equal(res1$eig,res2$eig)

## all(abs(res1$scores)-abs(res3$li)<1e-8)
## all(abs(res2$scores)-abs(res3$li)<1e-8)
## all(abs(res1$scores)-abs(res2$scores)<1e-8)

## all(abs(res1$loadings)-abs(res3$c1)<1e-8)
## all(abs(res2$loadings)-abs(res3$c1)<1e-8)
## all(abs(res1$loadings)-abs(res2$loadings)<1e-8)


## ## perform ordinary PCA
## titi <- dudi.pca(M, center=TRUE, scale=FALSE, scannf=FALSE, nf=4)


## ## check results
## all(round(abs(toto$scores), 10) == round(abs(titi$li), 10)) # MUST BE TRUE
## all.equal(toto$eig, titi$eig) # MUST BE TRUE
## all(round(abs(toto$loadings), 10)==round(abs(titi$c1), 10)) # MUST BE TRUE


## TEST WITH NAS ##
## M <- matrix(sample(c(0,1, NA), 1e5, replace=TRUE, prob=c(.495,.495,.01)), nrow=100)
## rownames(M) <- paste("ind", 1:100)

## ## perform glPca
## x <- new("genlight",M)
## toto <- glPca(x, nf=4)

## round(cor(toto$scores),10) # must be diag(1,4)
## round(t(toto$loadings) %*% toto$loadings,10) # must be diag(1,4)





## ## SPEED TESTS ##
## ## perform glPca
## M <- matrix(sample(c(0,1), 100*1e5, replace=TRUE), nrow=100)
## x <- new("genlight",M)
## system.time(titi <- dudi.pca(M,center=TRUE,scale=FALSE, scannf=FALSE, nf=4)) # 92 sec
## system.time(toto <- glPca(x, ,center=TRUE,scale=FALSE, useC=TRUE, nf=4)) # 102 sec

## M <- matrix(sample(c(0,1), 200*1e5, replace=TRUE), nrow=200)
## x <- new("genlight",M)
## system.time(titi <- dudi.pca(M,center=TRUE,scale=FALSE, scannf=FALSE, nf=4)) #  109 sec
## system.time(toto <- glPca(x, ,center=TRUE,scale=FALSE, useC=TRUE, nf=4)) #  360 sec

## M <- matrix(sample(c(0,1), 100*5e5, replace=TRUE), nrow=500)
## x <- new("genlight",M)
## system.time(titi <- dudi.pca(M,center=TRUE,scale=FALSE, scannf=FALSE, nf=4)) #  MEM LIMIT ISSUE
## system.time(toto <- glPca(x, ,center=TRUE,scale=FALSE, useC=TRUE, nf=4)) #  sec

## USE R PROFILING ##
## for glPca
## M <- matrix(sample(c(0,1), 100*1e5, replace=TRUE), nrow=100)
## x <- new("genlight",M)
## Rprof("glPca-prof.log")
## toto <- glPca(x, ,center=TRUE,scale=FALSE, useC=TRUE, nf=4) # 102 sec
## Rprof(NULL)
## res <- summaryRprof("glPca-prof.log")
## t <- res$by.total$total.time
## names(t) <- rownames(res$by.total)
## par(mar=c(7,4,4,2))
## barplot(t,las=3, cex.names=.7)


## ## for dudi.pca
## M <- matrix(sample(c(0,1), 100*1e5, replace=TRUE), nrow=100)
## Rprof("dudipca-prof.log")
## toto <- dudi.pca(M ,center=TRUE,scale=FALSE, scannf=FALSE, nf=4) # 102 sec
## Rprof(NULL)
## res <- summaryRprof("dudipca-prof.log")
## t <- res$by.total$total.time
## names(t) <- rownames(res$by.total)
## par(mar=c(7,4,4,2))
## barplot(t,las=3, cex.names=.7)


## test GLsum:
## library(adegenet)
## x <- new("genlight", lapply(1:50, function(i) sample(c(0,1,NA), 1e5, prob=c(.5, .49, .01), replace=TRUE)))
## res1 <- glSum(x, useC=FALSE)
## res2 <- glSum(x, useC=TRUE)
## res3 <- apply(as.matrix(x),2,sum,na.rm=TRUE)
## all(res1==res3) # must be TRUE
## all(res2==res3) # must be TRUE

## library(adegenet)
## x <- new("genlight", lapply(1:50, function(i) sample(c(0,1,2,NA), 1e5, prob=c(.5, .40, .09, .01), replace=TRUE)))
## res1 <- glSum(x, alleleAsUnit=FALSE, useC=FALSE)
## res2 <- glSum(x, alleleAsUnit=FALSE, useC=TRUE)
## res3 <- apply(as.matrix(x)/ploidy(x),2,sum,na.rm=TRUE)
## all.equal(res1,res3)
## all.equal(res2,res3)


## TEST PARALLELE C COMPUTATIONS IN GLDOTPROD
## system.time(toto <- glDotProd(x,multi=TRUE)) # 58 sec: cool!
## system.time(titi <- glDotProd(x,multi=FALSE)) # 245 sec
## all.equal(toto,titi)


## TEST PARALLELE C COMPUTATIONS IN GLPCA ##
## first dataset
## x <- new("genlight", lapply(1:50, function(i) sample(c(0,1,2,NA), 1e5, prob=c(.5, .40, .09, .01), replace=TRUE)))
## system.time(pca1 <- glPca(x, multi=FALSE, useC=FALSE, nf=1)) # no C, no parallel: 43 sec
## system.time(pca2 <- glPca(x, multi=FALSE, useC=TRUE, nf=1)) # just C: 248 sec
## system.time(pca3 <- glPca(x, multi=TRUE, useC=FALSE, nf=1, n.core=7)) # just parallel: 16 sec
## system.time(pca4 <- glPca(x, multi=TRUE, useC=TRUE, nf=1, n.core=7)) # C and parallel: 65 sec
## all.equal(pca1$scores^2, pca2$scores^2) # must be TRUE
## all.equal(pca1$scores^2, pca3$scores^2) # must be TRUE
## all.equal(pca1$scores^2, pca4$scores^2) # must be TRUE

## second dataset
## x <- new("genlight", lapply(1:500, function(i) sample(c(0,1,2,NA), 1e4, prob=c(.5, .40, .09, .01), replace=TRUE)))
## system.time(pca1 <- glPca(x, multi=FALSE, useC=FALSE, nf=1)) # no C, no parallel: 418 sec
## system.time(pca2 <- glPca(x, multi=FALSE, useC=TRUE, nf=1)) # just C:  496 sec
## system.time(pca3 <- glPca(x, multi=TRUE, useC=FALSE, nf=1, n.core=7)) # just parallel: 589 sec
## system.time(pca4 <- glPca(x, multi=TRUE, useC=TRUE, nf=1, n.core=7)) # C and parallel: 315 sec
## all.equal(pca1$scores^2, pca2$scores^2) # must be TRUE
## all.equal(pca1$scores^2, pca3$scores^2) # must be TRUE
## all.equal(pca1$scores^2, pca4$scores^2) # must be TRUE
