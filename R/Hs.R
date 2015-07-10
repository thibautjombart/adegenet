#'
#' Expected heterozygosity (Hs)
#'
#' This function computes the expected heterozygosity (Hs) within
#' populations of a \linkS4class{genpop} object. This function is
#' available for codominant markers (\code{@@type="codom"}) only. Hs is
#' commonly used for measuring within population genetic diversity (and
#' as such, it still has sense when computed from haploid data).
#'
#' @aliases Hs
#'
#' @rdname Hs
#'
#' @export
#'
#' @details
#' Let \emph{m(k)} be the number of alleles of locus \emph{k}, with a
#' total of \emph{K} loci. We note \eqn{f_i} the allele frequency of
#' allele \emph{i} in a given population. Then, \eqn{Hs} is given for a
#' given population by:\cr
#'
#' \eqn{\frac{1}{K} \sum_{k=1}^K (1 - \sum_{i=1}^{m(k)} f_i^2)}
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param x a \linkS4class{genind} or \linkS4class{genpop} object.
#' @param pop only used if x is a \linkS4class{genind}; an optional factor to be used as population; if not provided, pop(x) is used.
#'
#' @seealso \code{\link{Hs.test}} to test differences in Hs between two groups
#' @return a vector of Hs values (one value per population)
#' @examples
#' \dontrun{
#' data(nancycats)
#' Hs(genind2genpop(nancycats))
#' }
#'
Hs <- function(x, pop=NULL) {

    ## CHECKS
    if(is.genind(x)){
        if(!is.null(pop)) pop(x) <- pop
        x <- genind2genpop(x, quiet=TRUE)
    }
    if(!is.genpop(x)) stop("x is not a valid genpop object")
    if(x@type=="PA") stop("not implemented for presence/absence markers")


    ## MAIN COMPUTATIONS
    ## get number of typed loci per pop
    if(any(is.na(tab(x)))){
        temp <- apply(tab(x),1,tapply, locFac(x), function(e) !any(is.na(e)))
        nLoc.typed <- apply(temp,2,sum)
    } else {
        nLoc.typed <- nLoc(x)
    }

    F <- tab(x, freq=TRUE, NA.method="asis")
    res <- 1-(apply(F^2,1,sum, na.rm=TRUE))/nLoc.typed

    return(res)
} # end Hs


#'
#' Test differences in expected heterozygosity (Hs)
#'
#' This procedure permits to test if two groups have
#' significant differences in expected heterozygosity (Hs).
#' The test statistic used is simply the difference in Hs
#' between the two groups 'x' and 'y':
#'
#' \eqn{Hs(x) - Hs(y)}
#'
#' Individuals are randomly permuted between groups to obtain
#' a reference distribution of the test statistics.
#'
#' @rdname Hs.test
#' @aliases Hs.test
#'
#' @export
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param x a \linkS4class{genind} object.
#' @param y a \linkS4class{genind} object.
#' @param n.sim the number of permutations to be used to generate the reference distribution.
#' @param alter a character string indicating the alternative hypothesis
#'
#' @seealso \code{\link{Hs}} to compute Hs for different populations;
#' \code{\link[ade4]{as.randtest}} for the class of Monte Carlo tests.
#'
#' @return an object of the class randtest
#'
#' @examples
#' \dontrun{
#' data(microbov)
#' Hs(microbov)
#' test <- Hs.test(microbov[pop="Borgou"],
#'                 microbov[pop="Lagunaire"],
#'                 n.sim=499)
#' test
#' plot(test)
#' }
#'
Hs.test <- function(x, y, n.sim=999, alter=c("two-sided", "greater", "less")){
    ## CHECKS ##
    if(!is.genind(x)) stop("x is not a valid genind object")
    if(!is.genind(y)) stop("y is not a valid genind object")
    if(x@type=="PA" || y@type=="PA") stop("not implemented for presence/absence markers")
    alter <- match.arg(alter)

    ## POOL DATA ##
    xy <- repool(x,y)
    pop(xy) <- rep(c("x","y"), c(nInd(x), nInd(y)))

    ## AUXIL FUNCTION ##
    f1 <- function(x, pop) -diff(Hs(x, pop))

    ## COMPUTE STATS
    ## observed value
    obs <- f1(xy, pop(xy))

    ## permuted
    sim <- replicate(n.sim, f1(xy, sample(pop(xy))))

    ## randtest
    res <- as.randtest(sim=sim, obs=obs, alter=alter)
    res$call <- match.call()

    return(res)

} # end Hs.test
