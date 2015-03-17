

#' Compute scaled allele frequencies
#'
#' The generic function \code{scaleGen} is an analogue to the \code{scale}
#' function, but is designed with further arguments giving scaling options.\cr
#'
#' Methods are defined for \linkS4class{genind} and \linkS4class{genpop}
#' objects.  Both return data.frames of scaled allele frequencies.
#'
#'
#' @name scaleGen-methods
#' @aliases scaleGen scaleGen-methods scaleGen,genind-method
#' scaleGen,genpop-method
#' @docType methods
#' @param x a \linkS4class{genind} and \linkS4class{genpop} object
#' @param center a logical stating whether alleles frequencies should be
#' centred to mean zero (default to TRUE). Alternatively, a vector of numeric
#' values, one per allele, can be supplied: these values will be substracted
#' from the allele frequencies.
#' @param scale a logical stating whether alleles frequencies should be scaled
#' (default to TRUE). Alternatively, a vector of numeric values, one per
#' allele, can be supplied: these values will be substracted from the allele
#' frequencies.
#' @param truenames a logical indicating whether true labels (as opposed to
#' generic labels) should be used to name the output.
#' @param missing a character giving the treatment for missing values. Can be
#' "NA", "0" or "mean"
#' @return A matrix of scaled allele frequencies with genotypes
#' (\linkS4class{genind}) or populations in (\linkS4class{genpop}) in rows and
#' alleles in columns.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords methods manip
#' @examples
#'
#' \dontrun{
#' ## load data
#' data(microbov)
#' obj <- genind2genpop(microbov)
#'
#' ## compare different scaling
#' X1 <- scaleGen(obj)
#' X2 <- scaleGen(obj,met="bin")
#'
#' ## compute PCAs
#' pcaObj <- dudi.pca(obj,scale=FALSE,scannf=FALSE) # pca with no scaling
#' pcaX1 <- dudi.pca(X1,scale=FALSE,scannf=FALSE,nf=100) # pca with usual scaling
#' pcaX2 <- dudi.pca(X2,scale=FALSE,scannf=FALSE,nf=100) # pca with scaling for binomial variance
#'
#' ## get the loadings of alleles for the two scalings
#' U1 <- pcaX1$c1
#' U2 <- pcaX2$c1
#'
#'
#' ## find an optimal plane to compare loadings
#' ## use a procustean rotation of loadings tables
#' pro1 <- procuste(U1,U2,nf=2)
#'
#' ## graphics
#' par(mfrow=c(2,2))
#' # eigenvalues
#' barplot(pcaObj$eig,main="Eigenvalues\n no scaling")
#' barplot(pcaX1$eig,main="Eigenvalues\n usual scaling")
#' barplot(pcaX2$eig,main="Eigenvalues\n 'binomial' scaling")
#' # differences between loadings of alleles
#' s.match(pro1$scor1,pro1$scor2,clab=0,sub="usual -> binom (procustean rotation)")
#'
#' }
#'

setGeneric("scaleGen", function(x,...){standardGeneric("scaleGen")})

setMethod("scaleGen", "genind", function(x, center=TRUE, scale=TRUE,
                                         missing=c("NA","0","mean"), truenames=TRUE){

    THRES <- 1e-10
    missing <- match.arg(missing)
    ## checkType(x)

    ## handle "missing" arg
    if(missing %in% c("0","mean")){
        x <- na.replace(x, method=missing, quiet=TRUE)
    }

    ## handle specific cases
    if(scale[1]){
        ## get allele freq
        temp <- apply(x$tab,2,mean,na.rm=TRUE)
        if(x@type=="codom"){
            ## coerce sum of alleles freq to one (in case of missing data)
            temp <- tapply(temp, x$loc.fac, function(vec) return(vec/sum(vec)))
            pbar <- unlist(temp)
        } else {
            pbar <- temp
        }

        scale <- sqrt(pbar*(1-pbar))
    }

    X <- x$tab
    ## handle truenames
    if(truenames){
        X <- truenames(x)
        if(is.list(X)) { X <- X$tab }
    }

    ## return result
    res <- scale(X, center=center, scale=scale)

    ## issue a warning if some variances are null
    temp <- attr(res,"scaled:scale") < THRES
    if(any(temp)) {
        warning("Some scaling values are null.\n Corresponding alleles are removed.")
        res <- res[, !temp]
        attr(res,"scaled:center") <- attr(res,"scaled:center")[!temp]
        attr(res,"scaled:scale") <- attr(res,"scaled:scale")[!temp]
    }

    return(res)
})





setMethod("scaleGen", "genpop", function(x, center=TRUE, scale=TRUE,
                                         missing=c("NA","0","mean"), truenames=TRUE){

    THRES <- 1e-10
    missing <- match.arg(missing)
    ## checkType(x)

    ## make allele frequencies here
    if(x@type=="codom"){
        X <- makefreq(x,quiet=TRUE,missing=missing,truenames=truenames)$tab
    } else{
        X <- truenames(x) # keep binary data if type is PA
    }

    ## handle specific cases
    if(scale[1]){
        ## get allele freq
        temp <- apply(X,2,mean,na.rm=TRUE)
        if(x@type=="codom"){
            ## coerce sum of alleles freq to one (in case of missing data)
            temp <- tapply(temp, x$loc.fac, function(vec) return(vec/sum(vec)))
            pbar <- unlist(temp)
        } else{
            pbar <- temp
        }

        scale <- sqrt(pbar*(1-pbar))
    }

    ## return result

    res <- scale(X, center=center, scale=scale)

    ## issue a warning if some variances are null
    temp <- attr(res,"scaled:scale") < THRES
    if(any(temp)) {
        warning("Some scaling values are null.\n Corresponding alleles are removed.")
        res <- res[, !temp]
        attr(res,"scaled:center") <- attr(res,"scaled:center")[!temp]
        attr(res,"scaled:scale") <- attr(res,"scaled:scale")[!temp]
    }

    return(res)
})
