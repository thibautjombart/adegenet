#' Choose the number of clusters for genclust.em using BIC
#'
#' Do not use. We work on that stuff. Contact us if interested.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#'
#' @seealso \code{\link{genclust.em}} to generate individual clustering
#'     solutions, and \code{\link{BIC.genclust.em}} for computing BIC for
#'     \code{genclust.em} objects.
#'
#' @param max An integer indicating the maximum number of clusters to seek;
#'     \code{\link{genclust.em}} will be run for all k from 2 to max.
#'
#' @param IC A function computing the information criterion for
#'     \code{\link{BIC.genclust.em}} objects. Available statistics are
#'     \code{AIC} (default) and \code{BIC}.
#' 
#' @param IC.only A logical (TRUE by default) indicating if IC values only
#'     should be returned; if \code{FALSE}, full \code{genclust.em} objects are
#'     returned.
#'
#' @param ... Arguments passed to \code{\link{genclust.em}}.
#'
#' @seealso \code{\link{genclust.em} to generate clustering solutions,
#'     \code{\link{AIC.genclust.em}} for Akaike's criterion, and
#'     \code{\link{BIC.genclust.em}} for the Bayesian information criterion.
#' 
genclust.em.choose.k <- function(max, ..., IC = AIC, IC.only = TRUE) {

    ## This function is a glorified for loop which runs genclust.em for several
    ## values of 'k', from 2 to 'max'. It returns information criterion (AIC or
    ## BIC), and can also return the full genclust.em objects if needed. For
    ## k=1, AIC and BIC are computed via an internal (i.e. non-exported)
    ## procedure.
    
    max <- as.integer(max)
    if (any(!is.finite(max))) {
        stop("Values of k need to be finite.")
    }
    if (max < 2) {
        stop("maximum number of clusters should be at least 2")
    }
    k.values <- 2:max

    call.args <- list(...)
    genind.posi <- match("genind", sapply(call.args, class))
    if (is.na(genind.posi)) {
        stop("No genind provided in '...'.")
    }
    names(call.args)[genind.posi] <- "x"


    out.IC <- double(length(k.values))
    out.objects <- list(length(k.values))
    
    for (i in seq_along(k.values)) {
        ## get clustering solution for 'k'
        call.args$k <- k.values[i]
        out.objects[[i]] <- do.call(genclust.em, call.args)
    }

    names(out.objects) <- k.values
    out.IC <- .compute.null.IC(call.args$x)
    out.IC <- c(out.IC, vapply(out.objects, IC, double(1)))
    names(out.IC) <- 1:max
                 
    if (IC.only) {
        out <- out.IC
    } else {
        out <- list(AIC = out.IC, objects = out.objects)
        if (IC == BIC) {
            names(out)[1] <- "BIC"
        }
    }

    return(out)
}





## Non-exported procedure to compute the BIC for k = 1

## - 'x' is a genind object.
## - 'IC' is either AIC or BIC

.compute.null.IC <- function(x, IC = AIC) {
    group <- rep(1L, nInd(x))
    n.loc <- nLoc(x)
    genotypes <- tab(x)
    pop.freq <- tab(genind2genpop(x, pop = group, quiet = TRUE),
                    freq = TRUE)

    ## browser()

    ## get likelihoods of genotypes
    ll.mat <- apply(genotypes, 1, .ll.genotype, pop.freq, n.loc)
    ll.mat <- matrix(ll.mat, nrow = 1)
    
    ll <- .global.ll(group, ll.mat)

    ## make a fake genclust.em object to get BIC
    fake <- list(ll = ll,
                 group = group,
                 n.param = ncol(tab(x)) - nLoc(x))
    class(fake) <- "genclust.em"
    IC(fake)
}
