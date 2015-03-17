###########################
#
# Auxiliary functions for
# adegenet objects
#
# T. Jombart
###########################

##############################
# Method truenames for genind
##############################
setGeneric("truenames", function(x) standardGeneric("truenames"))

setMethod("truenames", signature(x="genind"), function(x){
    checkType(x)
    X <- x@tab
    if(!all(x@ind.names=="")) {rownames(X) <- x@ind.names}

    labcol <- locNames(x, withAlleles=TRUE)
    colnames(X) <- labcol

    if(!is.null(x@pop)){
        pop <- x@pop
        levels(pop) <- x@pop.names
        return(list(tab=X,pop=pop))
    }

    return(X)
}
          )





##############################
# Method truenames for genpop
##############################
setMethod("truenames",signature(x="genpop"), function(x){
    checkType(x)

    X <- x@tab
    if(!all(x@pop.names=="")) {rownames(X) <- x@pop.names}

    labcol <- locNames(x, withAlleles=TRUE)
    colnames(X) <- labcol

    return(X)
})







#' Separate data per locus
#'
#' The function \code{seploc} splits an object (\linkS4class{genind},
#' \linkS4class{genpop} or \linkS4class{genlight}) by marker. For
#' \linkS4class{genind} and \linkS4class{genpop} objects, the method returns a
#' list of objects whose components each correspond to a marker. For
#' \linkS4class{genlight} objects, the methods returns blocks of SNPs.
#'
#'
#' @name seploc
#' @aliases seploc seploc-methods seploc,ANY-method seploc,genind-method
#' seploc,genpop-method seploc,genlight-method
#' @docType methods
#' @param x a \linkS4class{genind} or a \linkS4class{genpop} object.
#' @param truenames a logical indicating whether true names should be used
#' (TRUE, default) instead of generic labels (FALSE).
#' @param res.type a character indicating the type of returned results, a
#' genind or genpop object (default) or a matrix of data corresponding to the
#' 'tab' slot.
#' @param n.block an integer indicating the number of blocks of SNPs to be
#' returned.
#' @param block.size an integer indicating the size (in number of SNPs) of the
#' blocks to be returned.
#' @param random should blocks be formed of contiguous SNPs, or should they be
#' made or randomly chosen SNPs.
#' @param parallel a logical indicating whether multiple cores -if available-
#' should be used for the computations (TRUE, default), or not (FALSE);
#' requires the package \code{parallel} to be installed.
#' @param n.cores if \code{parallel} is TRUE, the number of cores to be used in
#' the computations; if NULL, then the maximum number of cores available on the
#' computer is used.
#' @return The function \code{seploc} returns an list of objects of the same
#' class as the initial object, or a list of matrices similar to x\$tab.\cr
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{seppop}}, \code{\link{repool}}
#' @keywords manip
#' @examples
#'
#' \dontrun{
#' ## example on genind objects
#' data(microbov)
#'
#' # separate all markers
#' obj <- seploc(microbov)
#' names(obj)
#'
#' obj$INRA5
#'
#'
#' ## example on genlight objects
#' x <- glSim(100, 1000, 0, ploidy=2) # simulate data
#' x <- x[,order(glSum(x))] # reorder loci by frequency of 2nd allele
#' glPlot(x, main="All data") # plot data
#' foo <- seploc(x, n.block=3) # form 3 blocks
#' foo
#' glPlot(foo[[1]], main="1st block") # plot 1st block
#' glPlot(foo[[2]], main="2nd block") # plot 2nd block
#' glPlot(foo[[3]], main="3rd block") # plot 3rd block
#'
#' foo <- seploc(x, block.size=600, random=TRUE) # split data, randomize loci
#' foo # note the different block sizes
#' glPlot(foo[[1]])
#' }
#'

setGeneric("seploc", function(x, ...) standardGeneric("seploc"))

setMethod("seploc", signature(x="genind"), function(x,truenames=TRUE,res.type=c("genind","matrix")){
    if(x@type=="PA"){
        msg <- paste("seploc is not implemented for presence/absence markers")
        cat("\n",msg,"\n")
        return(invisible())
    }


    if(!is.genind(x)) stop("x is not a valid genind object")
    res.type <- match.arg(res.type)
    if(res.type=="genind") { truenames <- TRUE }

    temp <- x@loc.fac
    nloc <- length(levels(temp))
    levels(temp) <- 1:nloc

    kX <- list()

    for(i in 1:nloc){
        kX[[i]] <- matrix(x@tab[,temp==i],ncol=x@loc.nall[i])

        if(!truenames){
            rownames(kX[[i]]) <- rownames(x@tab)
            colnames(kX[[i]]) <- paste(names(x@loc.names)[i],names(x@all.names[[i]]),sep=".")
        }else{
            rownames(kX[[i]]) <- x@ind.names
            colnames(kX[[i]]) <- paste(x@loc.names[i],x@all.names[[i]],sep=".")
        }
    }

    if(truenames) {
        names(kX) <- x@loc.names
    } else{
        names(kX) <- names(x@loc.names)
    }

    prevcall <- match.call()
    if(res.type=="genind"){
        ## ploidy bug fixed by Zhian N. Kamvar
        ##kX <- lapply(kX, genind, pop=x@pop, prevcall=prevcall)
        kX <- lapply(kX, genind, pop=x@pop, prevcall=prevcall, ploidy=x@ploidy, type=x@type)
        for(i in 1:length(kX)){
            kX[[i]]@other <- x@other
        }
    }

    return(kX)
})



###########################
# Method seploc for genpop
###########################
setMethod("seploc", signature(x="genpop"), function(x,truenames=TRUE,res.type=c("genpop","matrix")){
     if(x@type=="PA"){
         msg <- paste("seploc is not implemented for presence/absence markers")
         cat("\n",msg,"\n")
         return(invisible())
    }


    if(!is.genpop(x)) stop("x is not a valid genpop object")
    res.type <- match.arg(res.type)
    if(res.type=="genpop") { truenames <- TRUE }

    temp <- x@loc.fac
    nloc <- length(levels(temp))
    levels(temp) <- 1:nloc

    kX <- list()

    for(i in 1:nloc){
        kX[[i]] <- matrix(x@tab[,temp==i],ncol=x@loc.nall[i])

        if(!truenames){
            rownames(kX[[i]]) <- rownames(x@tab)
            colnames(kX[[i]]) <- paste(names(x@loc.names)[i],names(x@all.names[[i]]),sep=".")
        }else{
            rownames(kX[[i]]) <- x@pop.names
            colnames(kX[[i]]) <- paste(x@loc.names[i],x@all.names[[i]],sep=".")
        }
    }

    if(truenames) {
        names(kX) <- x@loc.names
    } else{
        names(kX) <- names(x@loc.names)
    }

    prevcall <- match.call()
    if(res.type=="genpop"){
        kX <- lapply(kX, genpop, prevcall=prevcall, ploidy=x@ploidy, type=x@type)
        for(i in 1:length(kX)){
            kX[[i]]@other <- x@other
        }
    }

    return(kX)
})





###############
# '$' operator
###############
setMethod("$","genind",function(x,name) {
    return(slot(x,name))
})


setMethod("$<-","genind",function(x,name,value) {
   slot(x,name,check=TRUE) <- value
  return(x)
})









#' Separate genotypes per population
#'
#' The function \code{seppop} splits a \linkS4class{genind} or a
#' \linkS4class{genlight} object by population, returning a list of objects
#' whose components each correspond to a population.\cr
#'
#' For \linkS4class{genind} objects, the output can either be a list of
#' \linkS4class{genind} (default), or a list of matrices corresponding to the
#' \code{@@tab} slot.
#'
#'
#' @name seppop
#' @aliases seppop seppop-methods seppop,ANY-method seppop,genind-method
#' seppop,genlight-method
#' @docType methods
#' @param x a \linkS4class{genind} object
#' @param pop a factor giving the population of each genotype in 'x'. If not
#' provided, seeked in x\$pop.
#' @param truenames a logical indicating whether true names should be used
#' (TRUE, default) instead of generic labels (FALSE); used if res.type is
#' "matrix".
#' @param res.type a character indicating the type of returned results, a list
#' of \linkS4class{genind} object (default) or a matrix of data corresponding
#' to the 'tab' slots.
#' @param drop a logical stating whether alleles that are no longer present in
#' a subset of data should be discarded (TRUE) or kept anyway (FALSE, default).
#' @param treatOther a logical stating whether elements of the \code{@@other}
#' slot should be treated as well (TRUE), or not (FALSE). See details in
#' accessor documentations (\code{\link{pop}}).
#' @param quiet a logical indicating whether warnings should be issued when
#' trying to subset components of the \code{@@other} slot (TRUE), or not (FALSE,
#' default).
#' @param \dots further arguments passed to the genlight constructor.
#' @return According to 'res.type': a list of \linkS4class{genind} object
#' (default) or a matrix of data corresponding to the 'tab' slots.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{seploc}}, \code{\link{repool}}
#' @keywords manip
#' @examples
#'
#' \dontrun{
#' data(microbov)
#'
#' obj <- seppop(microbov)
#' names(obj)
#'
#' obj$Salers
#'
#'
#' #### example for genlight objects ####
#' x <- new("genlight", list(a=rep(1,1e3),b=rep(0,1e3),c=rep(1, 1e3)))
#' x
#'
#' pop(x) # no population info
#' pop(x) <- c("pop1","pop2", "pop1") # set population memberships
#' pop(x)
#' seppop(x)
#' as.matrix(seppop(x)$pop1)[,1:20]
#' as.matrix(seppop(x)$pop2)[,1:20,drop=FALSE]
#' }
#'

setGeneric("seppop", function(x, ...) standardGeneric("seppop"))

## genind
setMethod("seppop", signature(x="genind"), function(x,pop=NULL,truenames=TRUE,res.type=c("genind","matrix"), drop=FALSE, treatOther=TRUE, quiet=TRUE){
    ## checkType(x)

    ## misc checks
    if(!is.genind(x)) stop("x is not a valid genind object")
    if(is.null(pop)) { # pop taken from @pop
        if(is.null(x@pop)) stop("pop not provided and x@pop is empty")
        pop <- x@pop
        levels(pop) <- x@pop.names
    } else{
        pop <- factor(pop)
    }


    res.type <- match.arg(res.type)
    if(res.type=="genind") { truenames <- TRUE }

    ## pop <- x@pop # comment to take pop arg into account

    ## make a list of genind objects
    kObj <- lapply(levels(pop), function(lev) x[pop==lev, , drop=drop, treatOther=treatOther, quiet=quiet])
    names(kObj) <- levels(pop)

    ## res is a list of genind
    if(res.type=="genind"){ return(kObj) }

    ## res is list of matrices
    if(truenames) {
        res <- lapply(kObj, function(obj) truenames(obj)$tab)
    } else{
        res <- lapply(kObj, function(obj) obj$tab)
    }

    return(res)
}) # end seppop









#' Replace missing values (NA) from an object
#'
#' The generic function \code{na.replace} replaces NA in an object by
#' appropriate values as defined by the argument \code{method}.\cr
#'
#' Methods are defined for \linkS4class{genind} and \linkS4class{genpop}
#' objects.
#'
#' The argument "method" have the following effects:\cr - "0": missing values
#' are set to "0". An entity (individual or population) that is not type on a
#' locus has zeros for all alleles of that locus.\cr
#'
#' - "mean": missing values are set to the mean of the concerned allele,
#' computed on all available observations (without distinction of
#' population).\cr
#'
#' - "chi2": if a population is not typed for a marker, the corresponding count
#' is set to that of a theoretical count in of a Chi-squared test. This is
#' obtained by the product of the sums of both margins divided by the total
#' number of alleles.
#'
#' @name na.replace-methods
#' @aliases na.replace na.replace-methods na.replace,genind-method
#' na.replace,genpop-method
#' @docType methods
#' @param x a \linkS4class{genind} and \linkS4class{genpop} object
#' @param method a character string: can be "0" or "mean" for
#' \linkS4class{genind} objects, and "0" or "chi2" for \linkS4class{genpop}
#' objects.
#' @param quiet logical stating whether a message should be printed
#' (TRUE,default) or not (FALSE).
#' @return A \linkS4class{genind} and \linkS4class{genpop} object without
#' missing values.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords methods manip
#' @examples
#'
#' \dontrun{
#' data(nancycats)
#'
#' obj1 <- genind2genpop(nancycats)
#' # note missing data in this summary
#' summary(obj1)
#'
#' # NA are all in pop 17 and marker fca45
#' which(is.na(obj1$tab),TRUE)
#' truenames(obj1)[17,]
#'
#' # replace missing values
#' obj2 <- na.replace(obj1,"chi2")
#' obj2$loc.names
#'
#' # missing values where replaced
#' truenames(obj2)[,obj2$loc.fac=="L4"]
#' }
#'

setGeneric("na.replace", function(x, ...) standardGeneric("na.replace"))

## genind method
setMethod("na.replace", signature(x="genind"), function(x,method, quiet=FALSE){
    ## checkType(x)

    ## preliminary stuff
    validObject(x)
    if(!any(is.na(x@tab))) {
        if(!quiet) cat("\n Replaced 0 missing values \n")
        return(x)
    }
    method <- tolower(method)
    method <- match.arg(method, c("0","mean"))

    res <- x

    if(method=="0"){
        res@tab[is.na(x@tab)] <- 0
    }

    if(method=="mean"){
        f1 <- function(vec){
            m <- mean(vec,na.rm=TRUE)
            vec[is.na(vec)] <- m
            return(vec)
        }

        res@tab <- apply(x@tab, 2, f1)
    }

    if(!quiet){
        Nna <- sum(is.na(x@tab))
        cat("\n Replaced",Nna,"missing values \n")
    }

    return(res)

})




## genpop method
setMethod("na.replace", signature(x="genpop"), function(x,method, quiet=FALSE){
    ## checkType(x)

    ## preliminary stuff
    validObject(x)
    if(!any(is.na(x@tab))) {
        if(!quiet) cat("\n Replaced 0 missing values \n")
        return(x)
    }

    method <- tolower(method)
    method <- match.arg(method, c("0","chi2"))

    res <- x

    if(method=="0"){
        res@tab[is.na(x@tab)] <- 0
    }

    if(method=="chi2"){
        ## compute theoretical counts
        ## (same as in a Chi-squared)
        X <- x@tab
        sumPop <- apply(X,1,sum,na.rm=TRUE)
        sumLoc <- apply(X,2,sum,na.rm=TRUE)
        X.theo <- sumPop %o% sumLoc / sum(X,na.rm=TRUE)

        X[is.na(X)] <- X.theo[is.na(X)]
        res@tab <- X
    }

    if(!quiet){
        Nna <- sum(is.na(x@tab))
        cat("\n Replaced",Nna,"missing values \n")
    }

    return(res)
})





##################
# Function repool
##################


#' Pool several genotypes into a single dataset
#'
#' The function \code{repool} allows to merge genotypes from different
#' \linkS4class{genind} objects into a single 'pool' (i.e. a new
#' \linkS4class{genind}).  The markers have to be the same for all objects to
#' be merged, but there is no constraint on alleles.\cr
#'
#' This function can be useful, for instance, when hybrids are created using
#' \code{\link{hybridize}}, to merge hybrids with their parent population for
#' further analyses. Note that \code{repool} can also reverse the action of
#' \code{\link{seppop}}.
#'
#'
#' @param \dots can be i) a list whose components are valid
#' \linkS4class{genind} objects or, ii) several valid \linkS4class{genind}
#' objects separated by commas.
#' @return A \linkS4class{genind} object.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{seploc}}, \code{\link{seppop}}
#' @keywords manip
#' @examples
#'
#' \dontrun{
#' ## use the cattle breeds dataset
#' data(microbov)
#' temp <- seppop(microbov)
#' names(temp)
#'
#' ## hybridize salers and zebu -- nasty cattle
#' zebler <- hybridize(temp$Salers, temp$Zebu, n=40)
#' zebler
#'
#' ## now merge zebler with other cattle breeds
#' nastyCattle <- repool(microbov, zebler)
#' nastyCattle
#' }
#'
#' @export repool
repool <- function(...){

    ## preliminary stuff
    x <- list(...)
    if(is.list(x[[1]])) x <- x[[1]] ## if ... is a list, keep this list for x
    if(!inherits(x,"list")) stop("x must be a list")
    if(!all(sapply(x,is.genind))) stop("x is does not contain only valid genind objects")
    temp <- sapply(x,function(e) e$loc.names)
    if(!all(table(temp)==length(x))) stop("markers are not the same for all objects")
    temp <- sapply(x,function(e) e$ploidy)
    if(length(unique(temp)) != as.integer(1)) stop("objects have different levels of ploidy")



    ## extract info
    listTab <- lapply(x,genind2df,usepop=FALSE)
    getPop <- function(obj){
        if(is.null(obj$pop)) return(factor(rep(NA,nrow(obj$tab))))
        pop <- obj$pop
        levels(pop) <- obj$pop.names
        return(pop)
    }

    ## handle pop
    listPop <- lapply(x, getPop)
    pop <- unlist(listPop, use.names=FALSE)
    pop <- factor(pop)

    ## handle genotypes
    markNames <- colnames(listTab[[1]])
    listTab <- lapply(listTab, function(tab) tab[,markNames,drop=FALSE]) # resorting of the tabs

    ## bind all tabs by rows
    tab <- listTab[[1]]
    for(i in 2:length(x)){
        tab <- rbind(tab,listTab[[i]])
    }

    res <- df2genind(tab, pop=pop, ploidy=x[[1]]@ploidy, type=x[[1]]@type)
    res$call <- match.call()

    return(res)
} # end repool





#' Select genotypes of well-represented populations
#'
#' The function \code{selPopSize} checks the sample size of each population in
#' a \linkS4class{genind} object and keeps only genotypes of populations having
#' a given minimum size.
#'
#'
#' @name selPopSize
#' @aliases selPopSize selPopSize-methods selPopSize,ANY-method
#' selPopSize,genind-method
#' @docType methods
#' @param x a \linkS4class{genind} object
#' @param pop a vector of characters or a factor giving the population of each
#' genotype in 'x'. If not provided, seeked from x\$pop.
#' @param nMin the minimum sample size for a population to be retained. Samples
#' sizes strictly less than \code{nMin} will be discarded, those equal to or
#' greater than \code{nMin} are kept.
#' @return A \linkS4class{genind} object.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{seploc}}, \code{\link{repool}}
#' @keywords manip
#' @examples
#'
#' \dontrun{
#' data(microbov)
#'
#' table(pop(microbov))
#' obj <- selPopSize(microbov, n=50)
#'
#' obj
#' table(pop(obj))
#' }
#'

setGeneric("selPopSize", function(x, ...) standardGeneric("selPopSize"))

## genind method ##
setMethod("selPopSize", signature(x="genind"), function(x,pop=NULL,nMin=10){

    ## misc checks
    ## checkType(x)
    if(!is.genind(x)) stop("x is not a valid genind object")
    if(is.null(pop)) { # pop taken from @pop
        if(is.null(x@pop)) stop("pop not provided and x@pop is empty")
        pop <- x@pop
        levels(pop) <- x@pop.names
    } else{
        pop <- factor(pop)
    }

    ## select retained individuals
    effPop <- table(pop)
    popOk <- names(effPop)[effPop >= nMin]
    toKeep <- pop %in% popOk

    ## build result
    res <- x[toKeep]
    pop(res) <- pop[toKeep]

    return(res)
}) # end selPopSize








#' Assess polymorphism in genind/genpop objects
#'
#' The simple function \code{isPoly} can be used to check which loci are
#' polymorphic, or alternatively to check which alleles give rise to
#' polymorphism.
#'
#'
#' @name isPoly-methods
#' @aliases isPoly isPoly-methods isPoly,genind-method isPoly,genpop-method
#' @docType methods
#' @param x a \linkS4class{genind} and \linkS4class{genpop} object
#' @param by a character being "locus" or "allele", indicating whether results
#' should indicate polymorphic loci ("locus"), or alleles giving rise to
#' polymorphism ("allele").
#' @param thres a numeric value giving the minimum frequency of an allele
#' giving rise to polymorphism (defaults to 0.01).
#' @return A vector of logicals.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords methods manip
#' @examples
#'
#' \dontrun{
#' data(nancycats)
#' isPoly(nancycats,by="loc", thres=0.1)
#' isPoly(nancycats[1:3],by="loc", thres=0.1)
#' genind2df(nancycats[1:3])
#' }
#'

setGeneric("isPoly", function(x, ...) standardGeneric("isPoly"))

## genind method ##
setMethod("isPoly", signature(x="genind"), function(x, by=c("locus","allele"), thres=1/100){

    ## misc checks
    ## checkType(x)
    if(!is.genind(x)) stop("x is not a valid genind object")
    by <- match.arg(by)


    ## main computations ##

    ## PA case ##
    if(x@type=="PA") {
        allNb <- apply(x@tab, 2, mean, na.rm=TRUE) # allele frequencies
        toKeep <- (allNb >= thres) & (allNb <= (1-thres))
        return(toKeep)
    }


    ## codom case ##
    allNb <- apply(x@tab, 2, sum, na.rm=TRUE) # allele absolute frequencies

    if(by=="locus"){
        f1 <- function(vec){
            if(sum(vec) < 1e-10) return(FALSE)
            vec <- vec/sum(vec, na.rm=TRUE)
            if(sum(vec >= thres) >= 2) return(TRUE)
            return(FALSE)
        }

        toKeep <- tapply(allNb, x@loc.fac, f1)
    } else { # i.e. if mode==allele
        toKeep <- (allNb >= thres)
    }

    return(toKeep)
}) # end isPoly





## ## genpop method ##
## setMethod("isPoly", signature(x="genpop"), function(x, by=c("locus","allele"), thres=1/100){

##     ## misc checks
##     checkType(x)
##     if(!is.genpop(x)) stop("x is not a valid genind object")
##     by <- match.arg(by)


##     ## main computations ##
##     ## ## PA case ##
## ##     if(x@type=="PA") {
## ##         allNb <- apply(x@tab, 2, mean, na.rm=TRUE) # allele frequencies
## ##         toKeep <- (allNb >= thres) & (allNb <= (1-thres))
## ##         return(toKeep)
## ##     }


##     ## codom case ##
##     allNb <- apply(x@tab, 2, sum, na.rm=TRUE) # alleles absolute frequencies

##     if(by=="locus"){
##         f1 <- function(vec){
##             if(sum(vec) < 1e-10) return(FALSE)
##             vec <- vec/sum(vec, na.rm=TRUE)
##             if(sum(vec >= thres) >= 2) return(TRUE)
##             return(FALSE)
##         }

##         toKeep <- tapply(allNb, x@loc.fac, f1)
##     } else { # i.e. if mode==allele
##         toKeep <- allNb >= thres
##     }

##     return(toKeep)
## }) # end isPoly





######################
## miscellanous utils
######################

#######
# nLoc
#######
setGeneric("nLoc", function(x,...){
    standardGeneric("nLoc")
})



setMethod("nLoc","genind", function(x,...){
    return(length(x@loc.names))
})



setMethod("nLoc","genpop", function(x,...){
    return(length(x@loc.names))
})





#######
# nInd
#######
setGeneric("nInd", function(x,...){
    standardGeneric("nInd")
})



setMethod("nInd","genind", function(x,...){
    return(nrow(x@tab))
})





######
# pop
######
setGeneric("pop", function(x) {
  standardGeneric("pop")
})



setGeneric("pop<-",
           function(x, value) {
               standardGeneric("pop<-")
           })



setMethod("pop","genind", function(x){
    if(is.null(x@pop)) return(NULL)
    res <- x@pop
    levels(res) <- x@pop.names
    return(res)
})



setReplaceMethod("pop", "genind", function(x, value) {
    if(is.null(value)){
        x@pop <- NULL
        x@pop.names <- NULL
        return(x)
    }

    if(length(value) != nrow(x$tab)) stop("wrong length for population factor")

    ## coerce to factor (put levels in their order of appearance)
    newPop <- as.character(value)
    newPop <- factor(newPop, levels=unique(newPop))

    ## generic labels
    newPop.lab <- .genlab("P",length(levels(newPop)) )

    ## construct output
    x$pop.names <- levels(newPop)
    levels(newPop) <- newPop.lab
    x$pop <- newPop

    return(x)
})





###########
# locNames
###########
setGeneric("locNames", function(x,...){
    standardGeneric("locNames")
})

setGeneric("locNames<-", function(x, value) {
    standardGeneric("locNames<-")
})


setMethod("locNames","genind", function(x, withAlleles=FALSE, ...){
    ## return simply locus names
    if(x@type=="PA" | !withAlleles) return(x@loc.names)

    ## return locus.allele
    res <- rep(x@loc.names, x@loc.nall)
    res <- paste(res,unlist(x@all.names),sep=".")
    return(res)
})


setReplaceMethod("locNames","genind",function(x,value) {
    value <- as.character(value)
    if(length(value) != nLoc(x)) stop("Vector length does no match number of loci")
    names(value) <- names(locNames(x))
    slot(x,"loc.names",check=TRUE) <- value
    return(x)
})


setMethod("locNames","genpop", function(x, withAlleles=FALSE, ...){
    ## return simply locus names
    if(x@type=="PA" | !withAlleles) return(x@loc.names)

    ## return locus.allele
    res <- rep(x@loc.names, x@loc.nall)
    res <- paste(res,unlist(x@all.names),sep=".")
    return(res)
})


setReplaceMethod("locNames","genpop",function(x,value) {
    value <- as.character(value)
    if(length(value) != nLoc(x)) stop("Vector length does no match number of loci")
    names(value) <- names(locNames(x))
    slot(x,"loc.names",check=TRUE) <- value
    return(x)
})


###########
# indNames
###########
setGeneric("indNames", function(x,...){
    standardGeneric("indNames")
})

setGeneric("indNames<-", function(x, value){
    standardGeneric("indNames<-")
})

setMethod("indNames","genind", function(x, ...){
    return(x@ind.names)
})


setReplaceMethod("indNames","genind",function(x,value) {
    value <- as.character(value)
    if(length(value) != nInd(x)) stop("Vector length does no match number of individuals")
    names(value) <- names(indNames(x))
    slot(x,"ind.names",check=TRUE) <- value
    return(x)
})





##########
# alleles
##########
setGeneric("alleles", function(x,...){
    standardGeneric("alleles")
})

setGeneric("alleles<-", function(x, value){
    standardGeneric("alleles<-")
})

setMethod("alleles","genind", function(x, ...){
    return(x@all.names)
})

setReplaceMethod("alleles","genind", function(x, value){
    if(!is.list(value)) stop("replacement value must be a list")
    if(length(value)!=nLoc(x)) stop("replacement list must be of length nLoc(x)")
    if(any(sapply(value, length) != x$loc.nall)) stop("number of replacement alleles do not match that of the object")
    x@all.names <- value
    return(x)
})


setMethod("alleles","genpop", function(x, ...){
    return(x@all.names)
})

setReplaceMethod("alleles","genpop", function(x, value){
    if(!is.list(value)) stop("replacement value must be a list")
    if(length(value)!=nLoc(x)) stop("replacement list must be of length nLoc(x)")
    if(any(sapply(value, length) != x$loc.nall)) stop("number of replacement alleles do not match that of the object")
    x@all.names <- value
    return(x)
})



##########
## ploidy
##########
setGeneric("ploidy", function(x,...){
    standardGeneric("ploidy")
})

setGeneric("ploidy<-", function(x, value){
    standardGeneric("ploidy<-")
})

setMethod("ploidy","genind", function(x,...){
    return(x@ploidy)
})


setReplaceMethod("ploidy","genind",function(x,value) {
    value <- as.integer(value)
    if(any(value)<1) stop("Negative or null values provided")
    if(any(is.na(value))) stop("NA values provided")
    if(length(value)>1) warning("Several ploidy numbers provided; using only the first integer")
    slot(x,"ploidy",check=TRUE) <- value[1]
    return(x)
})


setMethod("ploidy","genpop", function(x,...){
    return(x@ploidy)
})


setReplaceMethod("ploidy","genpop",function(x,value) {
    value <- as.integer(value)
    if(any(value)<1) stop("Negative or null values provided")
    if(any(is.na(value))) stop("NA values provided")
    if(length(value)>1) warning("Several ploidy numbers provided; using only the first integer")
    slot(x,"ploidy",check=TRUE) <- value[1]
    return(x)
})






##########
## other
#########
setGeneric("other", function(x,...){
    standardGeneric("other")
})

setGeneric("other<-", function(x, value){
    standardGeneric("other<-")
})

setMethod("other","genind", function(x,...){
    if(length(x@other)==0) return(NULL)
    return(x@other)
})


setReplaceMethod("other","genind",function(x,value) {
    if( !is.null(value) && (!is.list(value) | is.data.frame(value)) ) {
        value <- list(value)
    }
    slot(x,"other",check=TRUE) <- value
    return(x)
})


setMethod("other","genpop", function(x,...){
    if(length(x@other)==0) return(NULL)
    return(x@other)
})


setReplaceMethod("other","genpop",function(x,value) {
    if( !is.null(value) && (!is.list(value) | is.data.frame(value)) ) {
        value <- list(value)
    }
    slot(x,"other",check=TRUE) <- value
    return(x)
})


