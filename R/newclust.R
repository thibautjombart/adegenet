#' Maximum-likelihood genetic clustering
#'
#' Do not use. We work on that stuff. Contact us if interested.  Current approach relies on
#' iterative EM algo for LL maximization and random group movements.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com} and Marie-Pauline Beugin
#'
#' @export
#'
#' @param x a \linkS4class{genind} object
#' @param k the number of clusters to look for
#' @param ini.group an optional factor defining the initial cluster configuration
#' @param max.em.iter the maximum number of iteration of the EM algorithm
#'
#' @examples
#' data(sim2pop)
#' x <- sim2pop
#' pop(x) <- sample(pop(x))
#'
#' ## check clustering when initial group is true one
#' clust <- newclust(sim2pop, 2, pop.ini=pop(sim2pop))
#' mean(clust$group==as.integer(pop(sim2pop)))
newclust <- function(x, k, ini.group=NULL, max.shake.iter=10, max.em.iter=100) {
    ## This function uses the EM algorithm to find ML group assignment of a set of genotypes stored
    ## in a genind object into 'k' clusters. We need an initial cluster definition to start with. The rest of the algorithm consists of:

    ## i) compute the matrix of allele frequencies
    ## ii) compute the likelihood of each genotype for each group
    ## iii) assign genotypes to the group for which they have the highest likelihood
    ## iv) go back to i) until convergence


    ## Set initial conditions: if initial pop is NULL, we create a random group definition (each
    ## clusters have same probability)
    if (is.null(ini.group)) {
        ini.group <- sample(seq_len(k), nInd(x), replace=TRUE)
    }

    ## make sure k and ini.group are compatible
    ini.group <- factor(ini.group)
    if (length(levels(ini.group)) != k) {
        stop("ini.group does not have k clusters")
    }

    ## initialisation
    group <- as.integer(ini.group) # groups are stored as integers
    genotypes <- tab(x)
    n.loc <- nLoc(x)
    n.ind <- nInd(x)
    n.groups <- length(unique(group))
    counter <- 1L
    converged <- FALSE

    ## get table of allele frequencies (columns) by population (rows)
    pop.freq <- tab(genind2genpop(x, pop=group, quiet=TRUE), freq=TRUE)

    ## get likelihoods of genotypes in every pop
    ll <- apply(genotypes, 1, ll.genotype, pop.freq, n.loc)

    ## assign individuals to their groups
    group <- apply(ll, 2, which.max)

    ## compute total log-likelihood for this clustering solution
    sum.ll <- global.ll(group, ll) # loglike of model



    ## The approach consists in alternating partial group randomisation (currently 20% of
    ## individuals assigned to a random group) and EM of group memberships optimizing the
    ## log-likelihood

    for (i in seq_len(max.shake.iter)) {
        ## compute current log-likelihood
        current.ll <- sum.ll

        ## shake groups: generate random groups for 10% of the data; move at least one individual
        ## each time

        n.to.move <- max(1, round(n.ind/10))
        to.move <- sample(1:n.ind, n.to.move, replace=FALSE)
        new.group <- group
        new.group[to.move] <- sample(seq_len(n.groups), size=n.to.move, replace=TRUE)

        ## EM-algorithm
        while(!converged && counter<=max.em.iter) {
            ## get table of allele frequencies (columns) by population (rows)
            pop.freq <- tab(genind2genpop(x, pop=group, quiet=TRUE), freq=TRUE)

            ll <- ll.all.genotypes(genotypes, pop.freq, n.loc)

            ## check convergence
            converged <- all(new.group == previous.group)
            counter <- counter + 1L

        }
        ## compute total log-likelihood for this clustering solution
        sum.ll <- global.ll(new.group, ll) # loglike of model
    }

    ## Metropolis algorithm to accept/reject new group
    p.accept <- exp(sum.ll - current.ll)
    if (runif(1) < p.accept) { # accept
        group <- new.group
        current.ll <- sum.ll
    } # if reject, we simply start again


    ## shape output and return
    proba <- round(prop.table(t(exp(ll)), 1), 2) # group membership proba


    out <- list(group = group,
                proba = proba,
                ll = sum.ll,
                converged = converged,
                n.iter = counter)
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




## Same function as above, except it is applied to all individuals in the dataset.
ll.all.genotypes <- function(genotypes, pop.freq, n.loc){
    apply(genotypes, 1, ll.genotype, pop.freq, n.loc)
}



## Non-exported function computing the total log-likelihood of the model given a vector of group
## assignments and a table of ll of genotypes in each group

global.ll <- function(group, ll){
    sum(t(ll)[cbind(seq_along(group), as.integer(group))])
}




## Non-exported function running one step of EM to find groups with optimum log-likelihood
## x: a genind object
## group: a factor / integer defining groups
assign.group.em1 <- function(x, group) {
    ## get table of allele frequencies (columns) by population (rows)
    pop.freq <- tab(genind2genpop(x, pop=group, quiet=TRUE), freq=TRUE)

    ## EXPECTATION: get likelihoods of genotypes in every group
    ll <- apply(genotypes, 1, ll.genotype, pop.freq, n.loc)

    ## MAXIMIZATION: assign individuals to most likely cluster
    apply(ll, 2, which.max)
}
