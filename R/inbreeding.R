

#' Likelihood-based estimation of inbreeding
#'
#' The function \code{inbreeding} estimates the inbreeding coefficient of an
#' individuals (F) by computing its likelihood function. It can return either
#' the density of probability of F, or a sample of F values from this
#' distribution. This operation is performed for all the individuals of a
#' \linkS4class{genind} object. Any ploidy greater than 1 is acceptable.
#'
#' Let \eqn{F} denote the inbreeding coefficient, defined as the probability
#' for an individual to inherit two identical alleles from a single ancestor.
#'
#' Let \eqn{p_i} refer to the frequency of allele \eqn{i} in the population.
#' Let \eqn{h} be an variable which equates 1 if the individual is homozygote,
#' and 0 otherwise. For one locus, the probability of being homozygote is
#' computed as:
#'
#' \eqn{ F + (1-F) \sum_i p_i^2}
#'
#' The probability of being heterozygote is: \eqn{1 - (F + (1-F) \sum_i p_i^2)}
#'
#' The likelihood of a genotype is defined as the probability of being the
#' observed state (homozygote or heterozygote). In the case of multilocus
#' genotypes, log-likelihood are summed over the loci.
#'
#' @aliases inbreeding
#' @param x an object of class \linkS4class{genind}.
#' @param pop a factor giving the 'population' of each individual. If NULL, pop
#' is seeked from \code{pop(x)}. Note that the term population refers in fact
#' to any grouping of individuals'.
#' @param truenames a logical indicating whether true names should be used
#' (TRUE, default) instead of generic labels (FALSE); used if res.type is
#' "matrix".
#' @param res.type a character string matching "sample", "function", or
#' "estimate" specifying whether the output should be a function giving the
#' density of probability of F values ("function"), the maximum likelihood
#' estimate of F from this distribution ("estimate"), or a sample of F values
#' taken from this distribution ("sample", default).
#' @param N an integer indicating the size of the sample to be taken from the
#' distribution of F values.
#' @param M an integer indicating the number of different F values to be used
#' to generate the sample. Values larger than N are recommended to avoid poor
#' sampling of the distribution.
#' @return A named list with one component for each individual, each of which
#' is a function or a vector of sampled F values (see \code{res.type}
#' argument).
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}\cr Zhian N.
#' Kamvar\cr
#' @seealso \code{\link{Hs}}: computation of expected heterozygosity.
#' @examples
#'
#' \dontrun{
#' ## cattle breed microsatellite data
#' data(microbov)
#'
#' ## isolate Lagunaire breed
#' lagun <- seppop(microbov)$Lagunaire
#'
#' ## estimate inbreeding - return sample of F values
#' Fsamp <- inbreeding(lagun, N=30)
#'
#' ## plot the first 10 results
#' invisible(sapply(Fsamp[1:10], function(e) plot(density(e), xlab="F",
#' xlim=c(0,1), main="Density of the sampled F values")))
#'
#' ## compute means for all individuals
#' Fmean=sapply(Fsamp, mean)
#' hist(Fmean, col="orange", xlab="mean value of F",
#' main="Distribution of mean F across individuals")
#'
#' ## estimate inbreeding - return proba density functions
#' Fdens <- inbreeding(lagun, res.type="function")
#'
#' ## view function for the first individual
#' Fdens[[1]]
#'
#' ## plot the first 10 functions
#' invisible(sapply(Fdens[1:10], plot, ylab="Density",
#' main="Density of probability of F values"))
#'
#' ## estimate inbreeding - return maximum likelihood estimates
#' Fest <- inbreeding(lagun, res.type = "estimate")
#' mostInbred <- which.max(Fest)
#' plot(Fdens[[mostInbred]], ylab = "Density", xlab = "F",
#'      main = paste("Probability density of F values\nfor", names(mostInbred)))
#' abline(v = Fest[mostInbred], col = "red", lty = 2)
#' legend("topright", legend = "MLE", col = "red", lty = 2)
#'
#' ## note that estimates and average samples are likely to be different.
#' plot(Fest, ylab = "F", col = "blue",
#'      main = "comparison of MLE and average sample estimates of F")
#' points(Fmean, pch = 2, col = "red")
#' arrows(x0 = 1:length(Fest), y0 = Fest,
#'        y1 = Fmean, x1 = 1:length(Fest), length = 0.125)
#' legend("topleft", legend = c("estimate", "sample"), col = c("blue", "red"),
#'        pch = c(1, 2), title = "res.type")
#' }
#'

###############
## inbreeding
###############
inbreeding <- function(x, pop=NULL, truenames=TRUE, res.type=c("sample","function","estimate"), N=200, M=N*10){
    ## CHECKS ##
    if(!is.genind(x)) stop("x is not a valid genind object")
    checkType(x)
    res.type <- match.arg(res.type)[1]

    ## if(x$ploidy != 2) stop("this inbreeding coefficient is designed for diploid genotypes only")
    PLO <- ploidy(x)

    if(!is.null(pop)) pop(x) <- pop
    if(is.null(x@pop) && is.null(pop)) {
        pop(x) <- factor(rep(1, nrow(x@tab)))
    }

      ## COMPUTATIONS ##

    ## get allele frequencies and \sum p_i^2 by pop and loc ##
    ## (generalized to any ploidy) ##
    ## tabfreq2 <- (makefreq(x = genind2genpop(x, quiet = TRUE), quiet=TRUE, truenames=truenames)$tab) ^2
    tabfreq2 <- (makefreq(x = genind2genpop(x, quiet = TRUE), quiet=TRUE, truenames=truenames)$tab) ^ PLO
    sumpi2 <- t(apply(tabfreq2, 1, tapply, x$loc.fac, sum))

    ## function to check a 1-locus genotype for homozigosity
    ## returns 1 if homoz, 0 otherwise
    ## !!! NOTE : reverse the values returned by f1 to obtain a strange thing !!!
    f1 <- function(gen){
        if(any(is.na(gen))) return(NA)
        if(any(round(gen, 10)==1)) return(1)
        return(0)
    }

    ## get the table of binary hetero/homo data
    if(truenames) {
        X <- truenames(x)$tab
    } else
    X <- x$tab

    homotab <- t(apply(X, 1, tapply, x@loc.fac, f1))


    ## get pi2 for the appropriate pop
    if(truenames){
    popx <- pop(x)
    } else {
        popx <- x$pop
    }

    popx <- as.character(popx)
    tabpi2 <- sumpi2[popx, , drop=FALSE]



    ## function returning a likelihood function - multi-locus ##
    ## x is a vector of 1/0 for each locus, with 1=homoz
    LIK <- function(x, sumpi2){
        ## return( sum(log(   x*(F+(1-F)*sumpi2) + (1-x)*(1-(F+(1-F)*sumpi2)) ) ,na.rm=TRUE) ) # returns the log-likelihood
        myEnv <- new.env()
        assign("x", x, envir=myEnv)
        assign("sumpi2", sumpi2, envir=myEnv)
        res <- function(F) {
            ## cat("\nx used:\n") # debugging
            ## print(x)
            ## cat("\nsumpi2 used:\n") # debugging
            ## print(sumpi2)

            return(exp(sum(log(   x*(F+(1-F)*sumpi2) + (1-x)*(1-(F+(1-F)*sumpi2)) ) ,na.rm=TRUE)))
        }
        environment(res) <- myEnv

        return(res)
    }


    ## get likelihood functions for all individuals
    res <- lapply(1:nrow(homotab), function(i) LIK(homotab[i,], tabpi2[i,]) )
    res <- lapply(res, Vectorize)

    ## name output
    if(truenames) {
        names(res) <- x$ind.names
    } else {
        names(res) <- names(x$tab)
    }


    ## IF WE RETURN FUNCTIONS ##
    if(res.type=="function"){
          return(res)
    }

    if (res.type == "estimate"){
      opfun <- function(x, ...) optimize(x, ...)[[1]]
      funval <- numeric(1)
      res <- vapply(res, FUN = opfun , FUN.VALUE = funval, interval = c(0, 1),
                    maximum = TRUE, tol = .Machine$double.eps^0.75)
      return(res)
    }


    ## IF WE RETURN SAMPLES ##
    ## function to get one sample
    getSample <- function(f){ # f is a vectorized density function
        x <- runif(M, 0, 1)
        fx <- f(x)
        fx <- fx - min(fx)
        if(sum(fx)<1e-14) {
            warning("Likelihood uniformly zero likely reflecting precision issue\nreturning uniformly distributed sample")
            return(runif(N))
        }
        fx <- fx/sum(fx)
        return(sample(x, size=N, prob=fx))
    }

    res <- lapply(res, getSample)
    return(res)
} # end inbreeding
