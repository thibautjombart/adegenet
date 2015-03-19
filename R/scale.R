
#' Compute scaled allele frequencies
#'
#' The generic function \code{scaleGen} is an analogue to the \code{scale}
#' function, but is designed with further arguments giving scaling options.\cr
#'
#' Methods are defined for \linkS4class{genind} and \linkS4class{genpop}
#' objects.  Both return data.frames of scaled allele frequencies.
#'
#'
#' @rdname scaleGen
#' @aliases scaleGen scaleGen-methods scaleGen,genind-method
#' scaleGen,genpop-method
#' @docType methods
#' @export
#' @param x a \linkS4class{genind} and \linkS4class{genpop} object
#' @param center a logical stating whether alleles frequencies should be
#' centred to mean zero (default to TRUE). Alternatively, a vector of numeric
#' values, one per allele, can be supplied: these values will be substracted
#' from the allele frequencies.
#' @param scale a logical stating whether alleles frequencies should be scaled
#' (default to TRUE). Alternatively, a vector of numeric values, one per
#' allele, can be supplied: these values will be substracted from the allele
#' frequencies.
#' @param truenames no longer used; kept for backward compatibility
#' @param NA.method a method to replace NA; asis: leave NAs as is; mean: replace by the mean allele frequencies; zero: replace by zero
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

#' @rdname scaleGen
#' @export
setMethod("scaleGen", "genind", function(x, center=TRUE, scale=TRUE,
                                         NA.method=c("asis","mean","zero"), truenames=TRUE){

    THRES <- 1e-10

    ## get table of frequencies
    out <- tab(x, NA.method=NA.method, freq=TRUE, quiet=TRUE)

    ## scale output
    out <- scale(out, center=center, scale=scale)

    ## issue a warning if some variances are null
    temp <- attr(out,"scaled:scale") < THRES
    if(any(temp)) {
        warning("Some scaling values are null.\n Corresponding alleles are removed.")
        out <- out[, !temp]
        attr(out,"scaled:center") <- attr(out,"scaled:center")[!temp]
        attr(out,"scaled:scale") <- attr(out,"scaled:scale")[!temp]
    }

    return(out)
})





#' @rdname scaleGen
#' @export
setMethod("scaleGen", "genpop", function(x, center=TRUE, scale=TRUE,
                                         NA.method=c("asis","mean","zero"), truenames=TRUE){

    THRES <- 1e-10

    ## get table of frequencies
    out <- tab(x, NA.method=NA.method, freq=TRUE, quiet=TRUE)

    ## scale output
    out <- scale(out, center=center, scale=scale)
    ## issue a warning if some variances are null
    temp <- attr(out,"scaled:scale") < THRES
    if(any(temp)) {
        warning("Some scaling values are null.\n Corresponding alleles are removed.")
        out <- out[, !temp]
        attr(out,"scaled:center") <- attr(out,"scaled:center")[!temp]
        attr(out,"scaled:scale") <- attr(out,"scaled:scale")[!temp]
    }

    return(out)
})

