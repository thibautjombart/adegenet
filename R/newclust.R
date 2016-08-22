#' Maximum-likelihood genetic clustering using EM algorithm
#'
#' Do not use. We work on that stuff. Contact us if interested.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com} and Marie-Pauline Beugin
#'
#' @export
#'
#' @param x a \linkS4class{genind} object
#' @param k the number of clusters to look for
#' @param pop.ini an optional factor defining the initial cluster configuration
#' @param max.iter the maximum number of iteration of the EM algorithm
#' @param detailed a logical stating whether extra details should be incorporated into the output;
#' these include group membership probability, indication of convergence, and the number of
#' iterations used before convergence
#'
#' @examples
#' data(sim2pop)
#' x <- sim2pop
#' pop(x) <- sample(pop(x))
#'
#' ## try function using true clusters as initial state
#'  genclust.em(x, 2, pop.ini = pop(sim2pop), detailed = FALSE)
genclust.em <- function(x, k, pop.ini = NULL, max.iter = 100, detailed = TRUE) {
    ## This function uses the EM algorithm to find ML group assignment of a set of genotypes stored
    ## in a genind object into 'k' clusters. We need an initial cluster definition to start with. The rest of the algorithm consists of:

    ## i) compute the matrix of allele frequencies
    ## ii) compute the likelihood of each genotype for each group
    ## iii) assign genotypes to the group for which they have the highest likelihood
    ## iv) go back to i) until convergence


    ## Set initial conditions: if initial pop is NULL, we create a random group definition (each
    ## clusters have same probability)
    if (is.null(pop.ini)) {
        pop.ini <- sample(seq_len(k), nInd(x), replace=TRUE)
    }

    ## make sure k and pop.ini are compatible
    pop.ini <- factor(pop.ini)
    if (length(levels(pop.ini)) != k) {
        stop("pop.ini does not have k clusters")
    }

    ## initialisation
    group <- factor(as.integer(pop.ini)) # levels are 1:k
    genotypes <- tab(x)
    n.loc <- nLoc(x)
    counter <- 1L
    converged <- FALSE

    while(!converged && counter<=max.iter) {

        ## get table of allele frequencies (columns) by population (rows)
        pop.freq <- tab(genind2genpop(x, pop=group, quiet=TRUE), freq=TRUE)

        ## get likelihoods of genotypes in every pop
        ll <- apply(genotypes, 1, ll.genotype, pop.freq, n.loc)

        ## assign individuals to most likely cluster
        previous.group <- group
        group <- apply(ll, 2, which.max)

        ## check convergence
        converged <- all(group == previous.group)
        counter <- counter + 1L

    }

    ## shape output and return
    out <- list(group = group, ll = global.ll(group, ll))

    if (detailed) {
        out$proba <- round(prop.table(t(exp(ll)), 1), 2) # group membership proba
        out$converged <- converged
        out$n.iter <- counter
    }

    return(out)
}





## Non-exported function which computes the log-likelihood of a genotype in every population. For
## now only works for diploid individuals. 'x' is a vector of allele counts; 'pop.freq' is a matrix
## of group allele frequencies, with groups in rows and alleles in columns.

## TODO: extend this to various ploidy levels, possibly optimizing procedures for haploids.

ll.genotype <- function(x, pop.freq, n.loc){
    ## homozygote (diploid)
    ## p(AA) = f(A)^2 for each locus
    ll.homoz <- apply(pop.freq, 1, function(f) sum(log(f[x==2L])) * 2)

    ## heterozygote (diploid, expl with 2 loci)
    ## p(Aa)p(Bb) = 2^n.loc * f(A)f(a) f(B)f(b)
    ll.heteroz <- apply(pop.freq, 1, function(f) sum(log(f[x==1L])) + n.loc * log(2))

    return(ll.homoz + ll.heteroz)
}






## Non-exported function computing the total log-likelihood of the model given a vector of group
## assignments and a table of ll of genotypes in each group

global.ll <- function(group, ll){
    sum(t(ll)[cbind(seq_along(group), as.integer(group))])
}

