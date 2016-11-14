#' Maximum-likelihood genetic clustering using EM algorithm
#'
#' Do not use. We work on that stuff. Contact us if interested.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com} and Marie-Pauline
#' Beugin
#'
#' @export
#'
#' @rdname genclust.em
#'
#' @param x a \linkS4class{genind} object
#'
#' @param k the number of clusters to look for
#'
#' @param pop.ini parameter indicating how the initial group membership should
#' be found. If \code{NULL}, groups are chosen at random, and the algorithm will
#' be run \code{n.start times}. If "kmeans", then the function
#' \code{find.clusters} is used to define initial groups. Alternatively, a
#' factor defining the initial cluster configuration can be provided.
#'
#' @param max.iter the maximum number of iteration of the EM algorithm
#'
#' @param n.start the number of times the EM algorithm is run, each time with
#' different random starting conditions
#'
#' @param hybrids a logical indicating if hybrids should be modelled
#' explicitely; this is currently implemented for 2 groups only.
#'
#' @param dim.ini the number of PCA axes to retain in the dimension reduction
#' step for \code{\link{find.clusters}}, if this method is used to define
#' initial group memberships (see argument \code{pop.ini}).
#'
#' @param hybrid.coef a vector of hybridization coefficients, defining the
#' proportion of hybrid gene pool coming from the first parental population;
#' this is symmetrized around 0.5, so that e.g. c(0.25, 0.5) will be converted
#' to c(0.25, 0.5, 0.75)
#'
#' @param parent.lab a vector of 2 character strings used to label the two
#' parental populations; only used if hybrids are detected (see argument
#' \code{hybrids})
#'
#' @param ... further arguments passed on to \code{\link{find.clusters}}
#'
#' @return
#'
#' The function \code{genclust.em} returns a list with the following
#' components:
#' \itemize{
#'
#' \item \code{$group} a factor indicating the maximum-likelihood assignment of
#' individuals to groups; if identified, hybrids are labelled after
#' hybridization coefficients, e.g. 0.5_A - 0.5_B for F1, 0.75_A - 0.25_B for
#' backcross F1 / A, etc.
#'
#' \item \code{$ll}: the log-likelihood of the model
#'
#' \item \code{$proba}: a matrix of group membership probabilities, with
#' individuals in rows and groups in columns; each value correspond to the
#' probability that a given individual genotype was generated under a given
#' group, under Hardy-Weinberg hypotheses.
#'
#' \item \code{$converged} a logical indicating if the algorithm converged; if
#' FALSE, it is doubtful that the result is an actual Maximum Likelihood
#' estimate.
#'
#' \item \code{$n.iter} an integer indicating the number of iterations the EM
#' algorithm was run for.
#'
#' }
#'
#' @examples
#' \dontrun{
#' data(microbov)
#'
#' ## try function using k-means initialization
#' grp.ini <- find.clusters(microbov, n.clust=15, n.pca=150)
#'
#' ## run EM algo
#' res <- genclust.em(microbov, 15, pop.ini = grp.ini$grp)
#' names(res)
#' res$converged
#' res$n.iter
#'
#' ## plot result
#' compoplot(res)
#'
#' ## flag potential hybrids
#' to.flag <- apply(res$proba,1,max)<.9
#' compoplot(res, subset=to.flag, show.lab=TRUE,
#'                  posi="bottomleft", bg="white")
#'
#'
#' ## Simulate hybrids F1
#' zebu <- microbov[pop="Zebu"]
#' salers <- microbov[pop="Salers"]
#' hyb <- hybridize(zebu, salers, n=30)
#' x <- repool(zebu, salers, hyb)
#'
#' ## method without hybrids
#' res.no.hyb <- genclust.em(x, k=2, hybrids=FALSE)
#' compoplot(res.no.hyb, col.pal=spectral, n.col=2)
#'
#' ## method with hybrids
#' res.hyb <- genclust.em(x, k=2, hybrids=TRUE)
#' compoplot(res.hyb, col.pal =
#'           hybridpal(col.pal = spectral), n.col = 2)
#'
#'
#' ## Simulate hybrids backcross (F1 / parental)
#' f1.zebu <- hybridize(hyb, zebu, 20, pop = "f1.zebu")
#' f1.salers <- hybridize(hyb, salers, 25, pop = "f1.salers")
#' y <- repool(x, f1.zebu, f1.salers)
#'
#' ## method without hybrids
#' res2.no.hyb <- genclust.em(y, k = 2, hybrids = FALSE)
#' compoplot(res2.no.hyb, col.pal = hybridpal(), n.col = 2)
#'
#' ## method with hybrids F1 only
#' res2.hyb <- genclust.em(y, k = 2, hybrids = TRUE)
#' compoplot(res2.hyb, col.pal = hybridpal(), n.col = 2)
#'
#' ## method with back-cross
#' res2.back <- genclust.em(y, k = 2, hybrids = TRUE, hybrid.coef = c(.25,.5))
#'  compoplot(res2.hyb, col.pal = hybridpal(), n.col = 2)
#'
#' }

genclust.em <- function(x, k, pop.ini = "kmeans", max.iter = 100, n.start=10,
                        hybrids = FALSE, dim.ini = 100,
                        hybrid.coef = NULL, parent.lab = c('A', 'B'), ...) {
    ## This function uses the EM algorithm to find ML group assignment of a set
    ## of genotypes stored in a genind object into 'k' clusters. We need an
    ## initial cluster definition to start with. The rest of the algorithm
    ## consists of:

    ## i) compute the matrix of allele frequencies

    ## ii) compute the likelihood of each genotype for each group

    ## iii) assign genotypes to the group for which they have the highest
    ## likelihood

    ## iv) go back to i) until convergence


    ## Disable multiple starts if the initial condition is not random
    use.random.start <- is.null(pop.ini)
    if (!use.random.start) {
        n.start <- 1L
    }

    if (n.start < 1L) {
        stop(sprintf(
            "n.start is less than 1 (%d); using n.start=1", n.start))
    }

    if (hybrids && k > 2) {
        warning(sprintf(
            "forcing k=2 for hybrid mode (requested k is %d)", k))
        k <- 2
    }


    ## Handle hybrid coefficients; these values reflect the contribution of the
    ## first parental population to the allele frequencies of the hybrid
    ## group. For instance, a value of 0.75 indicates that 'a' contributes to
    ## 75%, and 'b' 25% of the allele frequencies of the hybrid - a typical
    ## backcross F1 / a.

    if (hybrids) {
        if (is.null(hybrid.coef)) {
            hybrid.coef <- 0.5
        }
        hybrid.coef <- .tidy.hybrid.coef(hybrid.coef)
    }


    ## Initialisation using 'find.clusters'
    if (!is.null(pop.ini) &&
        tolower(pop.ini)[1] %in% c("kmeans", "k-means", "find.clusters")
        ) {
        pop.ini <- find.clusters(x, n.clust = k, n.pca = dim.ini, ...)$grp
    }

    ## There is one run of the EM algo for each of the n.start random initial
    ## conditions.
    ll <- -Inf # this will be the total loglike

    for (i in seq_len(n.start)) {

        ## Set initial conditions: if initial pop is NULL, we create a random
        ## group definition (each clusters have same probability)
        if (use.random.start) {
            pop.ini <- sample(seq_len(k), nInd(x), replace=TRUE)
        }

        ## process initial population, store levels
        pop.ini <- factor(pop.ini)
        lev.ini <- levels(pop.ini)[1:k] # k+1 would be hybrids

        ## ensure 'pop.ini' matches 'k'
        if (! (length(levels(pop.ini)) %in% c(k, k + length(hybrid.coef))) ) {
            stop("pop.ini does not have k clusters")
        }

        ## initialisation
        group <- factor(as.integer(pop.ini)) # set levels to 1:k (or k+1)
        genotypes <- tab(x)
        n.loc <- nLoc(x)
        counter <- 0L
        converged <- FALSE


        ## This is the actual EM algorithm

        while(!converged && counter<=max.iter) {

            ## get table of allele frequencies (columns) by population (rows);
            ## these are stored as 'pop.freq'; note that it will include extra
            ## rows for different types of hybrids too.

            if (hybrids) {
                pop(x) <- group
                id.parents <- .find.parents(x)
                x.parents <- x[id.parents]
                pop.freq <- tab(genind2genpop(x.parents, quiet=TRUE),
                                freq=TRUE)
                pop.freq <- rbind(pop.freq, # parents
                                  .find.freq.hyb(pop.freq, hybrid.coef)) # hybrids
            } else {
                pop.freq <- tab(genind2genpop(x, pop=group, quiet=TRUE),
                                freq=TRUE)
            }

            ## ensures no allele frequency is exactly zero
            pop.freq <- .tidy.pop.freq(pop.freq, locFac(x))

            ## get likelihoods of genotypes in every pop
            ll.mat <- apply(genotypes, 1, .ll.genotype, pop.freq, n.loc)

            ## assign individuals to most likely cluster
            previous.group <- group
            group <- apply(ll.mat, 2, which.max)

            ## check convergence
            ## converged <- all(group == previous.group)
            old.ll <- .global.ll(previous.group, ll.mat)
            new.ll <- .global.ll(group, ll.mat)
            if (!is.finite(new.ll)) {
                ## stop(sprintf("log-likelihood at iteration %d is not finite (%f)",
                ##              counter, new.ll))
            }
            converged <- abs(old.ll - new.ll) < 1e-14
            counter <- counter + 1L

        }

        ## ## store the best run so far
        ## new.ll <- .global.ll(group, ll.mat)

        if (new.ll > ll || i == 1L) {
            ## store results
            ll <- new.ll
            out <- list(group = group, ll = ll)

            ## group membership probability
            out$proba <- prop.table(t(exp(ll.mat)), 1)
            out$converged <- converged
            out$n.iter <- counter
        }
    } # end of the for loop

    ## restore labels of groups
    out$group <- factor(out$group)
    if (hybrids) {
        if (!is.null(parent.lab)) {
            lev.ini <- parent.lab
        }
        hybrid.labels <- paste0(hybrid.coef, "_", lev.ini[1], "-",
                               1 - hybrid.coef, "_", lev.ini[2])
        lev.ini <- c(lev.ini, hybrid.labels)
    }
    levels(out$group) <- lev.ini
        colnames(out$proba) <- lev.ini
    class(out) <- c("genclust.em", "list")
    return(out)
}











## Non-exported function which computes the log-likelihood of a genotype in
## every population. For now only works for diploid individuals. 'x' is a vector
## of allele counts; 'pop.freq' is a matrix of group allele frequencies, with
## groups in rows and alleles in columns.

## TODO: extend this to various ploidy levels, possibly optimizing procedures
## for haploids.

.ll.genotype <- function(x, pop.freq, n.loc){
    ## homozygote (diploid)
    ## p(AA) = f(A)^2 for each locus
    ll.homoz.one.indiv <- function(f) {
        sum(log(f[x == 2L]), na.rm = TRUE) * 2
    }

    ll.homoz <- apply(pop.freq, 1, ll.homoz.one.indiv)

    ## heterozygote (diploid, expl with 2 loci)
    ## p(Aa)p(Bb) = 2^n.loc * f(A)f(a) f(B)f(b)
    ll.hetero.one.indiv <- function(f) {
        sum(log(f[x == 1L]), na.rm = TRUE) + n.loc * log(2)
    }

    ll.heteroz <- apply(pop.freq, 1, ll.hetero.one.indiv)

    return(ll.homoz + ll.heteroz)
}






## Non-exported function computing the total log-likelihood of the model given a vector of group
## assignments and a table of ll of genotypes in each group

.global.ll <- function(group, ll){
    sum(t(ll)[cbind(seq_along(group), as.integer(group))], na.rm=TRUE)
}





## Non-exported function making a tidy vector of weights for allele frequencies
## of parental populations. It ensures that given any input vector of weights
## 'w' defining the types of hybrids, the output has the following properties:

## - strictly on ]0,1[

## - symmetric around 0.5, e.g. c(.25, .5) gives c(.25, .5, .75)

## - sorted by decreasing values (i.e. hybrid types are sorted by decreasing
## proximity to the first parental population.

.tidy.hybrid.coef <- function(w) {
    w <- w[w > 0 & w < 1]
    w <- sort(unique(round(c(w, 1-w), 4)), decreasing = TRUE)
    w
}





## Non-exported function determining vectors of allele frequencies in hybrids
## from 2 parental populations. Different types of hybrids are determined by
## weights given to the allele frequencies of the parental populations. Only one
## such value is provided and taken to be the weight of the 1st parental
## population; the complementary frequency is derived for the second parental
## population.

## Parameters are:

## - x: matrix of allele frequencies for population 'a' (first row) and 'b'
## (second row), where allele are in columns.


## - w: a vector of weights for 'a' and 'b', each value determining a type of
## hybrid. For instance, 0.5 is for F1, 0.25 for backcrosses F1/parental, 0.125
## for 2nd backcross F1/parental, etc.

## The output is a matrix of allele frequencies with hybrid types in rows and
## alleles in columns.

.find.freq.hyb <- function(x, w) {
    out <- cbind(w, 1-w) %*% x
    rownames(out) <- w
    out
}




## Non-exported function trying to find the two parental populations in a genind
## object containing 'k' clusters. The parental populations are defined as the
## two most distant clusters. The other clusters are deemed to be various types
## of hybrids. The output is a vector of indices identifying the individuals
## from the parental populations.

.find.parents <- function(x) {
    ## matrix of pairwise distances between clusters, using Nei's distance
    D <- as.matrix(dist.genpop(genind2genpop(x, quiet = TRUE), method = 1))
    parents <- which(abs(max(D)-D) < 1e-14, TRUE)[1,]
    out <- which(as.integer(pop(x)) %in% parents)
    out
}





## Non-exported function enforcing a minimum allele frequency in a table of
## allele frequency. As we are not accounting for the uncertainty in allele
## frequencies, we need to allow for genotypes to be generated from a population
## which does not have the genotype's allele represented, even if this is at a
## low probability. The transformation is ad-hoc, and has the form:
##
## g(f_i) = (a + f_i /  \sum(a + f_i))

## where f_i is the i-th frequency in a given locus. However, this ensures that
## the output has two important properties:

## - it sums to 1
## - it contains no zero

## By default, we set 'a' to 0.01.

## Function inputs are:

## - 'pop.freq': matrix of allele frequencies, with groups in rows and alleles in
## columns

## - 'loc.fac': a factor indicating which alleles belong to which locus, as
## returned by 'locFac([a genind])'

.tidy.pop.freq <- function(pop.freq, loc.fac) {
    g <- function(f, a = .01) {
        (a + f) / sum(a + f)
    }

    out <- matrix(unlist(apply(pop.freq, 1, tapply, loc.fac, g),
                          use.names = FALSE),
                  byrow=TRUE, nrow=nrow(pop.freq))
    dimnames(out) <- dimnames(pop.freq)
    return(out)
}
