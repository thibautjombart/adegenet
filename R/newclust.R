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
#' @param n.start the number of times the EM algorithm is run, each time with different random
#' starting conditions
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
genclust.em <- function(x, k, pop.ini = NULL, max.iter = 100, n.start=10, detailed = TRUE) {
    ## This function uses the EM algorithm to find ML group assignment of a set of genotypes stored
    ## in a genind object into 'k' clusters. We need an initial cluster definition to start with. The rest of the algorithm consists of:

    ## i) compute the matrix of allele frequencies
    ## ii) compute the likelihood of each genotype for each group
    ## iii) assign genotypes to the group for which they have the highest likelihood
    ## iv) go back to i) until convergence


    ## Disable multiple starts if the initial condition is not random
    use.random.start <- is.null(pop.ini)
    if (!use.random.start) {
        n.start <- 1L
    }

    if (use.random.start && n.start < 1L) {
        warning(sprintf("n.start is less than 1 (%d); using n.start=1", n.start))
        n.start <- 1L
    }

    ## There is one run of the EM algo for each of the n.start random initial conditions.
    ll <- -Inf # this will be the total loglike

    for (i in seq_len(n.start)) {

        ## Set initial conditions: if initial pop is NULL, we create a random group definition (each
        ## clusters have same probability)
        if (use.random.start) {
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
            ll.mat <- apply(genotypes, 1, ll.genotype, pop.freq, n.loc)

            ## assign individuals to most likely cluster
            previous.group <- group
            group <- apply(ll.mat, 2, which.max)

            ## check convergence
            converged <- all(group == previous.group)
            counter <- counter + 1L

        }

        ## store the best run so far
        new.ll <- global.ll(group, ll.mat)

        if (new.ll > ll) {
            ## store results
            ll <- new.ll
            out <- list(group = group, ll = ll)

            if (detailed) {
                out$proba <- round(prop.table(t(exp(ll.mat)), 1), 2) # group membership proba
                out$converged <- converged
                out$n.iter <- counter
            }

        }
    } # end of the for loop


    return(out)
}





#' EM-MCMC algorithm for genetic clustering
#'
#' Do not use.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#' @export
#'
#' @rdname emmcmc
#'
#' @inheritParams genclust.em
#' @param n.iter the number of iterations of the MCMC
#' @param sample.every the frequency of sampling from the MCMC
#' @param max.em.iter the maximum number of iterations for the EM algorithm
#' @param prop.move the proportion of individuals moved across groups at each MCMC iteration
#'
#' @examples
#' data(sim2pop)
#' hyb <- hybridize(sim2pop[pop=1], sim2pop[pop=2], n=10)
#' x <- repool(sim2pop, hyb)
#'

## This algorithm relies on using EM within steps of a MCMC based on the likelihood (no prior) or
## observed genotypes given their group memberships and the group allele frequencies.
genclust.emmcmc <- function(x, k, n.iter = 100, sample.every = 10, pop.ini = NULL,
                            max.em.iter = 20, prop.move = 0.1, detailed) {

    ## initialize the algorithm
    mcmc <- vector(length = floor(n.iter / sample.every) + 1, mode="list")

    ## the object handled is a list with a $group (group membership) and a $ll (total loglike)
    mcmc[[1]] <- current.state <- genclust.em(x = x, k = k, pop.ini = pop.ini,
                             max.iter = max.em.iter, n.start = 10,
                             detailed = FALSE)
    mcmc[[1]]$step <- 1L
    counter <- 1L
    n.accept <- 0

    group.lev <- seq_len(k)

    ## proceed to the MCMC
    for (i in seq_len(n.iter)) {
        ## move groups
        new.groups <- move.groups(current.state$group, group.lev, prop.move = prop.move)

        ## compute loglike difference
        new.state <- genclust.em(x = x, k = k, pop.ini = new.groups,
                             max.iter = max.em.iter, detailed = FALSE)

        ## accept/reject (Metropolis algorithm)
        if(log(runif(1)) <= (new.state$ll - current.state$ll)) { ## accept
            current.state <- new.state
            n.accept <- n.accept + 1
        } # reject is implicitly: do nothing

        ## save result if needed
        if (i %% sample.every == 0) {
            counter <- counter + 1
            mcmc[[counter]] <- current.state
            mcmc[[counter]]$step <- i
        }
    } # end of MCMC

    ## shape result and return output
    out <- data.frame(step = sapply(mcmc, function(e) e$step),
                      ll = sapply(mcmc, function(e) e$ll))
    out.groups <- t(sapply(mcmc, function(e) e$group))
    out <- cbind.data.frame(out, out.groups)

    ## print(sprintf("\nProportion of movements accepted: %f (%d out of %d)", n.accept/n.iter, n.accept, n.iter))
    class(out) <- c("emmcmc", "data.frame")
    return(out)

}




#' @rdname emmcmc
#' @export
#'
plot.emmcmc <- function(x, y = "ll", type = c("trace", "hist", "density", "groups"), burnin = 0){
    ## CHECKS ##
    type <- match.arg(type)
    if (!y %in% names(x)) {
        stop(paste(y,"is not a column of x"))
    }

    ## GET DATA TO PLOT ##
    if (burnin > max(x$step)) {
        stop("burnin exceeds the number of steps in x")
    }
    x <- x[x$step>burnin,,drop=FALSE]

    ## MAKE PLOT ##
    if (type == "trace") {
        out <- ggplot(x) + geom_line(aes_string(x="step", y=y)) +
            labs(x="Iteration", y=y, title=paste("trace:",y))
    }
    if (type=="hist") {
        out <- ggplot(x) + geom_histogram(aes_string(x=y)) +
            geom_point(aes_string(x=y, y=0), shape="|", alpha=0.5, size=3) +
                labs(x=y, title=paste("histogram:",y))
    }

    if (type=="density") {
        out <- ggplot(x) + geom_density(aes_string(x=y)) +
            geom_point(aes_string(x=y, y=0), shape="|", alpha=0.5, size=3) +
                labs(x=y, title=paste("density:",y))
    }

    if (type=="groups") {
        groups <- x[, -c(1,2), drop=FALSE]
        out.dat <- data.frame(individual = as.vector(col(groups)),
                              group = factor(unlist(groups)))

        ## horrible kludge to standardize ylab
        out.expr <- paste("ggplot(out.dat, aes(x=individual))",
                          sprintf("geom_bar(aes(fill=group, y = (..count..)/%s))", nrow(x)),
                          "labs(y = 'assignement probability')",
                          sep=" + ")
        out <- eval(parse(text = out.expr))
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



## Non-exported function randomly moving group memberships for some individuals.
## -'x': a vector of group membership
## - 'pool': a vector of integer 1:K where K is the number of groups
## - 'prop.move': the proportion of individuals to move
move.groups <- function(x, pool, prop.move = 0.2) {
    n.to.move <- max(1, round(prop.move * length(x)))
    new.groups <- sample(pool, n.to.move, replace = TRUE)
    x[sample(seq_along(x), n.to.move, replace = FALSE)] <- new.groups
    return(x)
}
