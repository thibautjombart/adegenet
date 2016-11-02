## #' EM-MCMC algorithm for genetic clustering
## #'
## #' Do not use.
## #'
## #' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
## #' @export
## #'
## #' @rdname emmcmc
## #'
## #' @inheritParams genclust.em
## #' @param n.iter the number of iterations of the MCMC
## #' @param sample.every the frequency of sampling from the MCMC
## #' @param max.em.iter the maximum number of iterations for the EM algorithm
## #' @param prop.move the proportion of individuals moved across groups at each MCMC iteration
## #' @param hybrids a logical indicating if hybrids should be modelled explicitely; this is currently
## #' implemented for 2 groups only.
## #'
## #' @examples
## #' \dontrun{
## #' ## run analysis
## #' data(sim2pop)
## #' res <- genclust.emmcmc(sim2pop, k=2, n.iter=1e3, sample.every=50)
## #' plot(res)
## #'
## #' ## get summary
## #' smry <- summary(res, burn=250)
## #' compoplot(t(smry$proba))
## #' table(smry$group, pop(sim2pop))
## #'
## #' ## other example, generating hybrids
## #' p1 <- hybridize(sim2pop[pop=1], sim2pop[pop=1], n=50, pop="p1")
## #' p2 <- hybridize(sim2pop[pop=2], sim2pop[pop=2], n=50, pop="p2")
## #' hyb <- hybridize(p1, p2, n=50)
## #' x <- repool(p1, p2, hyb)
## #'
## #' ## make a pca
## #' pca1 <- dudi.pca(tab(x), scale=FALSE, scannf=FALSE)
## #' s.class(pca1$li, pop(x))
## #'
## #' ## run MCMC
## #' res <- genclust.emmcmc(x, k=2, n.iter=1e3, sample.every=50, hybrids=TRUE)
## #' plot(res)
## #'
## #' ## get summary
## #' smry <- summary(res, burn=250)
## #' compoplot(t(smry$proba), col.pal=spectral,
## #' n.col=2, txt.leg=paste("group", 1:2))
## #'
## #' table(smry$group, pop(x))
## #'
## #' ## same analysis, initialized with k-means
## #' clust.ini <- find.clusters(x, n.clust=2, n.pca=20)$grp
## #' res <- genclust.emmcmc(x, k=2, pop.ini=clust.ini)
## #' plot(res) # it has converged
## #' plot(res, type="group") # show groups
## #' }

## ## This algorithm relies on using EM within steps of a MCMC based on the
## ## likelihood (no prior) or observed genotypes given their group memberships and
## ## the group allele frequencies.

## genclust.emmcmc <- function(x, k, n.iter = 100, sample.every = 10,
##                             pop.ini = NULL, max.em.iter = 20, prop.move = 0.1,
##                             hybrids = FALSE, detailed = FALSE) {

##     ## initialize the algorithm
##     mcmc <- vector(length = floor(n.iter / sample.every) + 1, mode="list")

##     ## the object handled is a list with a $group (group membership) and a $ll
##     ## (total loglike)
##     mcmc[[1]] <- genclust.em(x = x, k = k, pop.ini = pop.ini,
##                              max.iter = max.em.iter, n.start = 10,
##                              hybrids = hybrids, detailed = FALSE)

##     current.state <- mcmc[[1]]

##     mcmc[[1]]$step <- 1L
##     counter <- 1L
##     n.accept <- 0

##     group.lev <- seq_len(k)

##     ## proceed to the MCMC
##     for (i in seq_len(n.iter)) {
##         ## move groups
##         new.groups <- .move.groups(current.state$group, prop.move = prop.move)

##         ## compute loglike difference
##         new.state <- genclust.em(x = x, k = k, pop.ini = new.groups,
##                                  max.iter = max.em.iter,
##                                  hybrids = hybrids, detailed = FALSE)

##         ## accept/reject (Metropolis algorithm)
##         if(log(runif(1)) <= (new.state$ll - current.state$ll)) { ## accept
##             current.state <- new.state
##             n.accept <- n.accept + 1
##         } # reject is implicitly: do nothing

##         ## save result if needed
##         if (i %% sample.every == 0) {
##             counter <- counter + 1
##             mcmc[[counter]] <- current.state
##             mcmc[[counter]]$step <- i
##         }
##     } # end of MCMC

##     ## shape result and return output
##     out <- data.frame(step = sapply(mcmc, function(e) e$step),
##                       ll = sapply(mcmc, function(e) e$ll))
##     out.groups <- t(sapply(mcmc, function(e) e$group))
##     out <- cbind.data.frame(out, out.groups)

##     ## print(sprintf("\nProportion of movements accepted: %f (%d out of %d)",
##     ## n.accept/n.iter, n.accept, n.iter))
##     class(out) <- c("genclust.emmcmc", "data.frame")
##     return(out)

## }



## #' @rdname emmcmc
## #' @export
## #' @param object a 'emmcmc' object
## summary.genclust.emmcmc <- function(object, burnin = 0, ...) {
##     if (burnin > max(object$step)) {
##         stop(sprintf(
##             "Burnin (%d) exceeds chain length (%d)", burnin, max(object$step)))
##     }
##     groups <- object[object$step > burnin , -(1:2)]
##     n.lev <- length(unique(unlist(groups)))
##     out <- list()
##     out$proba <- apply(groups, 2, function(e)
##                        vapply(seq_len(n.lev),
##                               function(i) mean(e==i), double(1)))
##     out$group <- apply(out$proba, 2, which.max)
##     return(out)
## }



## #' @rdname emmcmc
## #' @export
## #'
## #' @param y the data to plot
## #'
## #' @param type a character string indicating the type of graph to produce: the
## #' trace ('trace'), a histogram ('hist'), a density plot ('density') or a plot
## #' of group composition ('groups')
## #'
## #' @param burnin the number of iterations of the MCMC to discard
## #'
## #' @param ... further arguments passed to other methods
## #'
## plot.genclust.emmcmc <- function(x, y = "ll",
##                                  type = c("trace", "hist", "density", "groups"),
##                                  burnin = 0, ...){
##     ## CHECKS ##
##     type <- match.arg(type)
##     if (!y %in% names(x)) {
##         stop(paste(y,"is not a column of x"))
##     }

##     ## GET DATA TO PLOT ##
##     if (burnin > max(x$step)) {
##         stop("burnin exceeds the number of steps in x")
##     }
##     x <- x[x$step>burnin,,drop=FALSE]

##     ## MAKE PLOT ##
##     if (type == "trace") {
##         out <- ggplot2::ggplot(x) +
##                ggplot2::geom_line(ggplot2::aes_string(x="step", y=y)) +
##                ggplot2::labs(x="Iteration", y=y, title=paste("trace:",y))
##     }
##     if (type=="hist") {
##         out <- ggplot2::ggplot(x) +
##             ggplot2::geom_histogram(ggplot2::aes_string(x=y)) +
##             ggplot2::geom_point(ggplot2::aes_string(x=y, y=0),
##                                 shape="|", alpha=0.5, size=3) +
##                 ggplot2::labs(x=y, title=paste("histogram:",y))
##     }

##     if (type=="density") {
##         out <- ggplot2::ggplot(x) +
##             ggplot2::geom_density(ggplot2::aes_string(x=y)) +
##             ggplot2::geom_point(ggplot2::aes_string(x=y, y=0),
##                                 shape="|", alpha=0.5, size=3) +
##                 ggplot2::labs(x=y, title=paste("density:",y))
##     }

##     if (type=="groups") {
##         groups <- x[, -c(1,2), drop=FALSE]
##         out.dat <- data.frame(individual = as.vector(col(groups)),
##                               group = factor(unlist(groups)))

##         ## horrible kludge to standardize ylab
##         bar.txt <- paste0("ggplot2::geom_bar(",
##                           "ggplot2::aes(fill=group, y = (..count..)/%s))")
##         out.expr <- paste(
##             "ggplot2::ggplot(out.dat, aes(x=individual))",
##             sprintf(bar.txt, nrow(x)),
##             "labs(y = 'assignement probability')", sep=" + ")

##         out <- eval(parse(text = out.expr)) }

##     return(out)
## }






## ## Non-exported function randomly moving group memberships for some individuals.
## ## -'x': a factor indicating group membership
## ## - 'prop.move': the proportion of individuals to move

## .move.groups <- function(x, prop.move = 0.2) {
##     n.to.move <- max(1, round(prop.move * length(x)))
##     id.to.move <- sample(seq_along(x), n.to.move, replace = FALSE)
##     new.groups <- sample(levels(x), n.to.move, replace = TRUE)
##     x[id.to.move] <- new.groups
##     return(x)
## }
