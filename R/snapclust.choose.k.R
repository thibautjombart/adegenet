#' Choose the number of clusters for snapclust using AIC, BIC or AICc
#'
#' This function implements methods for investigating the optimal number of
#' genetic clusters ('k') using the fast maximum-likelihood genetic clustering
#' approach described in Beugin et al (2018). The method runs
#' \link{\code{snapclust}} for varying values of 'k', and computes the requested
#' summary statistics for each clustering solution to assess goodness of
#' fit. The method is fully documented in a dedicated tutorial which can be
#' accessed using \code{adegenetTutorial("snapclust")}.
#' 
#' @details The method is described in Beugin et al (2018) A fast likelihood
#'   solution to the genetic clustering problem. Methods in Ecology and
#'   Evolution \url{https://doi.org/10.1111/2041-210X.12968}. A dedicated
#'   tutorial is available by typing \code{adegenetTutorial("snapclust")}.
#' 
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#'
#' @seealso \code{\link{snapclust}} to generate individual clustering solutions,
#' and \code{\link{BIC.snapclust}} for computing BIC for \code{snapclust}
#' objects.
#'
#' @param max An integer indicating the maximum number of clusters to seek;
#'     \code{\link{snapclust}} will be run for all k from 2 to max.
#'
#' @param IC A function computing the information criterion for
#'     \code{\link{snapclust}} objects. Available statistics are
#'     \code{AIC} (default), \code{AICc}, and \code{BIC}.
#'
#' @param IC.only A logical (TRUE by default) indicating if IC values only
#'     should be returned; if \code{FALSE}, full \code{snapclust} objects are
#'     returned.
#'
#' @param ... Arguments passed to \code{\link{snapclust}}.
#'
snapclust.choose.k <- function(max, ..., IC = AIC, IC.only = TRUE) {

    ## This function is a glorified for loop which runs snapclust for several
    ## values of 'k', from 2 to 'max'. It returns information criterion (AIC or
    ## BIC), and can also return the full snapclust objects if needed. For
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
    genind.posi <- which(vapply(call.args, inherits, logical(1), "genind"))
    if (length(genind.posi) == 0) {
        stop("No genind provided in '...'.")
    } else if (length(genind.posi) > 1){
        warning("Too many genind objects provided in '...'. I am taking the first one.")
        genind.posi <- genind.posi[1]
    }
    names(call.args)[genind.posi] <- "x"


    out.IC <- double(length(k.values))
    out.objects <- list(length(k.values))

    for (i in seq_along(k.values)) {
        ## get clustering solution for 'k'
        call.args$k <- k.values[i]
        out.objects[[i]] <- do.call(snapclust, call.args)
    }

    names(out.objects) <- k.values
    out.IC <- .compute.null.IC(call.args$x)
    out.IC <- c(out.IC, vapply(out.objects, IC, double(1)))
    names(out.IC) <- 1:max

    if (IC.only) {
        out <- out.IC
    } else {
        out <- list(AIC = out.IC, objects = out.objects)

        ## names stat as appropriate
        names(out)[1] <- deparse(substitute(IC))
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

    if (!is.genind(x)) {
      stop("x is not a valid genind object")
    }

    if (any(ploidy(x) > 2)) {
      stop("snapclust not currently implemented for ploidy > 2")
    }

    if (all(ploidy(x) == 1)) {
      .ll.genotype <- .ll.genotype.haploid
    } else if (all(ploidy(x) == 2)) {
      .ll.genotype <- .ll.genotype.diploid
    } else {
      stop("snapclust not currently implemented for varying ploidy")
    }

    ## get likelihoods of genotypes
    ll.mat <- apply(genotypes, 1, .ll.genotype, pop.freq, n.loc)
    ll.mat <- matrix(ll.mat, nrow = 1)

    ll <- .global.ll(group, ll.mat)

    ## make a fake snapclust object to get IC
    fake <- list(ll = ll,
                 group = group,
                 n.param = ncol(tab(x)) - nLoc(x))
    class(fake) <- "snapclust"
    IC(fake)
}
