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
#'
newclust <- function(x, k) {
    ## This function uses the EM algorithm to find ML group assignment of a set of genotypes stored
    ## in a genind object into 'k' clusters. We need an initial cluster definition to start with. The rest of the algorithm consists of:

    ## i) compute the matrix of allele frequencies
    ## ii) compute the likelihood of each genotype for each group
    ## iii) assign genotypes to the group for which they have the highest likelihood
    ## iv) go back to i) until convergence
}


