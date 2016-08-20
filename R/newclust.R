#' Maximum-likelihood genetic clustering
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
#'
newclust <- function(x, k, pop.ini=NULL) {
    ## This function uses the EM algorithm to find ML group assignment of a set of genotypes stored
    ## in a genind object into 'k' clusters. We need an initial cluster definition to start with. The rest of the algorithm consists of:

    ## i) compute the matrix of allele frequencies
    ## ii) compute the likelihood of each genotype for each group
    ## iii) assign genotypes to the group for which they have the highest likelihood
    ## iv) go back to i) until convergence


    ## Set initial conditions: if initial pop is NULL, we use the pop available in 'x'; if 'x' has
    ## no pop, we create a random group definition (each clusters have same probability)
    if (is.null(pop.ini)) {
        if(!is.null(pop(x))) {
            pop.ini <- pop(x)
        }
    } else {
        pop.ini <- sample(seq_len(k), nInd(x), replace=TRUE)
    }

    ## make sure k and pop.ini are compatible
    pop.ini <- factor(pop.ini)
    if (length(levels(pop.ini)) != k) {
        stop("pop.ini does not have k clusters")
    }

    group <- pop.ini
    genotypes <- tab(x)
    n.loc <- nLoc(x)

    ## while(!converged) {
    ## get table of allele frequencies (columns) by population (rows)
    pop.freq <- tab(genind2genpop(x, pop=group, quiet=TRUE), freq=TRUE)

    ## get likelihoods of genotypes in every pop
    ll <- apply(genotypes, 1, ll.genotype, pop.freq, n.loc)
    ## }
}





## Non-exported function which computes the log-likelihood of a genotype in every population. For
## now only works for diploid individuals. 'x' is a vector of allele counts; 'pop.freq' is a matrix
## of group allele frequencies, with groups in rows and alleles in columns.

## TODO: extend this to various ploidy levels, possibly optimizing procedures for haploids.
ll.genotype <- function(x, pop.freq, n.loc){
    pop.freq <- t(pop.freq)

    ## homozygote (diploid)
    ## p(AA) = f(A)^2 for each locus
    ll.homoz <- apply(pop.freq, 1, function(f) sum(log(f[x==2L])) * 2)

    ## heterozygote (diploid, expl with 2 loci)
    ## p(Aa)p(Bb) = 2^n.loc * f(A)f(a) f(B)f(b)
    ll.heteroz <- apply(pop.freq, 1, function(f) sum(log(f[x==1L])) + n.loc * log(2))

    return(ll.homoz + ll.heteroz)
}


