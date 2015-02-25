#######
## dapc
########


#' Discriminant Analysis of Principal Components (DAPC)
#' 
#' These functions implement the Discriminant Analysis of Principal Components
#' (DAPC, Jombart et al. 2010). This method descibes the diversity between
#' pre-defined groups. When groups are unknown, use \code{find.clusters} to
#' infer genetic clusters. See 'details' section for a succint description of
#' the method, and \code{vignette("adegenet-dapc")} for a tutorial. Graphical
#' methods for DAPC are documented in \code{\link{scatter.dapc}} (see
#' \code{?scatter.dapc}).
#' 
#' \code{dapc} is a generic function performing the DAPC on the following types
#' of objects:\cr - \code{data.frame} (only numeric data)\cr - \code{matrix}
#' (only numeric data)\cr - \code{\linkS4class{genind}} objects (genetic
#' markers)\cr - \code{\linkS4class{genlight}} objects (genome-wide SNPs)
#' 
#' These methods all return an object with class \code{dapc}.
#' 
#' Functions that can be applied to these objects are (the ".dapc" can be
#' ommitted):
#' 
#' - \code{print.dapc}: prints the content of a \code{dapc} object.\cr -
#' \code{summary.dapc}: extracts useful information from a \code{dapc}
#' object.\cr - \code{predict.dapc}: predicts group memberships based on DAPC
#' results.\cr - \code{xvalDapc}: performs cross-validation of DAPC using
#' varying numbers of PCs (and keeping the number of discriminant functions
#' fixed); it currently has methods for \code{data.frame} and \code{matrix}.\cr
#' 
#' DAPC implementation calls upon \code{\link[ade4]{dudi.pca}} from the
#' \code{ade4} package (except for \linkS4class{genlight} objects) and
#' \code{\link[MASS]{lda}} from the \code{MASS} package. The \code{predict}
#' procedure uses \code{\link[MASS]{predict.lda}} from the \code{MASS} package.
#' 
#' \code{as.lda} is a generic with a method for \code{dapc} object which
#' converts these objects into outputs similar to that of \code{lda.default}.
#' 
#' The Discriminant Analysis of Principal Components (DAPC) is designed to
#' investigate the genetic structure of biological populations. This
#' multivariate method consists in a two-steps procedure. First, genetic data
#' are transformed (centred, possibly scaled) and submitted to a Principal
#' Component Analysis (PCA). Second, principal components of PCA are submitted
#' to a Linear Discriminant Analysis (LDA). A trivial matrix operation allows
#' to express discriminant functions as linear combination of alleles,
#' therefore allowing one to compute allele contributions. More details about
#' the computation of DAPC are to be found in the indicated reference.
#' 
#' DAPC does not infer genetic clusters ex nihilo; for this, see the
#' \code{\link{find.clusters}} function.
#' 
#' @aliases dapc dapc.data.frame dapc.matrix dapc.genind dapc.dudi
#' dapc.genlight print.dapc summary.dapc predict.dapc as.lda as.lda.dapc
#' @param x \code{a data.frame}, \code{matrix}, or \code{\linkS4class{genind}}
#' object. For the \code{data.frame} and \code{matrix} arguments, only
#' quantitative variables should be provided.
#' @param grp,pop a \code{factor} indicating the group membership of
#' individuals; for \code{scatter}, an optional grouping of individuals.
#' @param n.pca an \code{integer} indicating the number of axes retained in the
#' Principal Component Analysis (PCA) step. If \code{NULL}, interactive
#' selection is triggered.
#' @param n.da an \code{integer} indicating the number of axes retained in the
#' Discriminant Analysis step. If \code{NULL}, interactive selection is
#' triggered.
#' @param center a \code{logical} indicating whether variables should be
#' centred to mean 0 (TRUE, default) or not (FALSE). Always TRUE for
#' \linkS4class{genind} objects.
#' @param scale a \code{logical} indicating whether variables should be scaled
#' (TRUE) or not (FALSE, default). Scaling consists in dividing variables by
#' their (estimated) standard deviation to account for trivial differences in
#' variances.
#' @param var.contrib a \code{logical} indicating whether the contribution of
#' original variables (alleles, for \linkS4class{genind} objects) should be
#' provided (TRUE, default) or not (FALSE). Such output can be useful, but can
#' also create huge matrices when there is a lot of variables.
#' @param pca.info a \code{logical} indicating whether information about the
#' prior PCA should be stored (TRUE, default) or not (FALSE). This information
#' is required to predict group membership of new individuals using
#' \code{predict}, but makes the object slightly bigger.
#' @param pca.select a \code{character} indicating the mode of selection of PCA
#' axes, matching either "nbEig" or "percVar". For "nbEig", the user has to
#' specify the number of axes retained (interactively, or via \code{n.pca}).
#' For "percVar", the user has to specify the minimum amount of the total
#' variance to be preserved by the retained axes, expressed as a percentage
#' (interactively, or via \code{perc.pca}).
#' @param perc.pca a \code{numeric} value between 0 and 100 indicating the
#' minimal percentage of the total variance of the data to be expressed by the
#' retained axes of PCA.
#' @param list() further arguments to be passed to other functions. For
#' \code{dapc.matrix}, arguments are to match those of \code{dapc.data.frame};
#' for \code{dapc.genlight}, arguments passed to \code{\link{glPca}}
#' @param glPca an optional \code{\link{glPca}} object; if provided, dimension
#' reduction is not performed (saving computational time) but taken directly
#' from this object.
#' @param object a \code{dapc} object.
#' @param truenames a \code{logical} indicating whether true (i.e.,
#' user-specified) labels should be used in object outputs (TRUE, default) or
#' not (FALSE).
#' @param dudi optionally, a multivariate analysis with the class \code{dudi}
#' (from the ade4 package). If provided, prior PCA will be ignored, and this
#' object will be used as a prior step for variable orthogonalisation.
#' @param newdata an optional dataset of individuals whose membership is
#' seeked; can be a data.frame, a matrix, a \linkS4class{genind} or a
#' \linkS4class{genlight} object, but object class must match the original
#' ('training') data. In particular, variables must be exactly the same as in
#' the original data. For \linkS4class{genind} objects, see
#' \code{\link{repool}} to ensure matching of alleles.
#' @param prior,dimen,method see \code{?predict.lda}.
#' @return === dapc objects ===\cr The class \code{dapc} is a list with the
#' following components:\cr \item{call}{the matched call.} \item{n.pca}{number
#' of PCA axes retained} \item{n.da}{number of DA axes retained}
#' \item{var}{proportion of variance conserved by PCA principal components}
#' \item{eig}{a numeric vector of eigenvalues.} \item{grp}{a factor giving
#' prior group assignment} \item{prior}{a numeric vector giving prior group
#' probabilities} \item{assign}{a factor giving posterior group assignment}
#' \item{tab}{matrix of retained principal components of PCA}
#' \item{loadings}{principal axes of DAPC, giving coefficients of the linear
#' combination of retained PCA axes.} \item{ind.coord}{principal components of
#' DAPC, giving the coordinates of individuals onto principal axes of DAPC;
#' also called the discriminant functions.} \item{grp.coord}{coordinates of the
#' groups onto the principal axes of DAPC.} \item{posterior}{a data.frame
#' giving posterior membership probabilities for all individuals and all
#' clusters.} \item{var.contr}{(optional) a data.frame giving the contributions
#' of original variables (alleles in the case of genetic data) to the principal
#' components of DAPC.} \item{match.prp}{a list, where each item is the
#' proportion of individuals correctly matched to their original population in
#' cross-validation.}
#' 
#' === other outputs ===\cr Other functions have different outputs:\cr -
#' \code{summary.dapc} returns a list with 6 components: \code{n.dim} (number
#' of retained DAPC axes), \code{n.pop} (number of groups/populations),
#' \code{assign.prop} (proportion of overall correct assignment),
#' \code{assign.per.pop} (proportion of correct assignment per group),
#' \code{prior.grp.size} (prior group sizes), and \code{post.grp.size}
#' (posterior group sizes), \code{xval.dapc}, \code{xval.genind} and
#' \code{xval} (all return a list of four lists, each one with as many items as
#' cross-validation runs.  The first item is a list of \code{assign}
#' components, the secon is a list of \code{posterior} components, the thirs is
#' a list of \code{ind.score} components and the fourth is a list of
#' \code{match.prp} items, i.e. the prortion of the validation set correctly
#' matched to its original population)
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \itemize{ \item \code{\link{xvalDapc}}: selection of the optimal
#' numbers of PCA axes retained in DAPC using cross-validation.
#' 
#' \item \code{\link{scatter.dapc}}, \code{\link{assignplot}},
#' \code{\link{compoplot}}: graphics for DAPC.
#' 
#' \item \code{\link{find.clusters}}: to identify clusters without prior.
#' 
#' \item \code{\link{dapcIllus}}: a set of simulated data illustrating the DAPC
#' 
#' \item \code{\link{eHGDP}}, \code{\link{H3N2}}: empirical datasets
#' illustrating DAPC }
#' @references Jombart T, Devillard S and Balloux F (2010) Discriminant
#' analysis of principal components: a new method for the analysis of
#' genetically structured populations. BMC Genetics11:94.
#' doi:10.1186/1471-2156-11-94
#' @keywords multivariate
#' @examples
#' 
#' ## data(dapcIllus), data(eHGDP), and data(H3N2) illustrate the dapc
#' ## see ?dapcIllus, ?eHGDP, ?H3N2
#' ##
#' \dontrun{
#' example(dapcIllus)
#' example(eHGDP)
#' example(H3N2)
#' }
#' 
#' ## H3N2 EXAMPLE ##
#' data(H3N2)
#' pop(H3N2) <- factor(H3N2$other$epid)
#' dapc1 <- dapc(H3N2, var.contrib=FALSE, scale=FALSE, n.pca=150, n.da=5)
#' 
#' ## remove internal segments and ellipses, different pch, add MStree
#' scatter(dapc1, cell=0, pch=18:23, cstar=0, mstree=TRUE, lwd=2, lty=2)
#' 
#' ## only ellipse, custom labels
#' scatter(dapc1, cell=2, pch="", cstar=0, posi.da="top",
#' lab=paste("year\n",2001:2006), axesel=FALSE, col=terrain.colors(10))
#' 
#' 
#' ## SHOW COMPOPLOT ON MICROBOV DATA ##
#' data(microbov)
#' dapc1 <- dapc(microbov, n.pca=20, n.da=15)
#' compoplot(dapc1, lab="")
#' 
#' 
#' 
#' 
#' \dontrun{
#' ## EXAMPLE USING GENLIGHT OBJECTS ##
#' ## simulate data
#' x <- glSim(50,4e3-50, 50, ploidy=2)
#' x
#' plot(x)
#' 
#' ## perform DAPC
#' dapc1 <- dapc(x, n.pca=10, n.da=1)
#' dapc1
#' 
#' ## plot results
#' scatter(dapc1, scree.da=FALSE)
#' 
#' ## SNP contributions
#' loadingplot(dapc1$var.contr)
#' loadingplot(tail(dapc1$var.contr, 100), main="Loading plot - last 100 SNPs")
#' 
#' 
#' 
#' ## USE "PREDICT" TO PREDICT GROUPS OF NEW INDIVIDUALS ##
#' ## load data
#' data(sim2pop)
#' 
#' ## we make a dataset of:
#' ## 30 individuals from pop A
#' ## 30 individuals from pop B
#' ## 30 hybrids
#' 
#' ## separate populations and make F1
#' temp <- seppop(sim2pop)
#' temp <- lapply(temp, function(e) hybridize(e,e,n=30)) # force equal popsizes
#' 
#' ## make hybrids
#' hyb <- hybridize(temp[[1]], temp[[2]], n=30)
#' 
#' ## repool data - needed to ensure allele matching
#' newdat <- repool(temp[[1]], temp[[2]], hyb)
#' pop(newdat) <- rep(c("pop A", "popB", "hyb AB"), c(30,30,30))
#' 
#' ## perform the DAPC on the first 2 pop (60 first indiv)
#' dapc1 <- dapc(newdat[1:60],n.pca=5,n.da=1)
#' 
#' ## plot results
#' scatter(dapc1, scree.da=FALSE)
#' 
#' ## make prediction for the 30 hybrids
#' hyb.pred <- predict(dapc1, newdat[61:90])
#' hyb.pred
#' 
#' ## plot the inferred coordinates (circles are hybrids)
#' points(hyb.pred$ind.scores, rep(.1, 30))
#' 
#' ## look at assignment using assignplot
#' assignplot(dapc1, new.pred=hyb.pred)
#' title("30 indiv popA, 30 indiv pop B, 30 hybrids")
#' 
#' ## image using compoplot
#' compoplot(dapc1, new.pred=hyb.pred, ncol=2)
#' title("30 indiv popA, 30 indiv pop B, 30 hybrids")
#' 
#' ## CROSS-VALIDATION ##
#' data(sim2pop)
#' xval <- xvalDapc(sim2pop@tab, pop(sim2pop), n.pca.max=100, n.rep=3)
#' xval
#' boxplot(xval$success~xval$n.pca, xlab="Number of PCA components",
#' ylab="Classification succes", main="DAPC - cross-validation")
#' 
#' }
#' 
#' 
#' 
#' @export dapc
dapc <- function (x, ...) UseMethod("dapc")

###################
## dapc.data.frame
###################
dapc.data.frame <- function(x, grp, n.pca=NULL, n.da=NULL,
                            center=TRUE, scale=FALSE, var.contrib=TRUE, pca.info=TRUE,
                            pca.select=c("nbEig","percVar"), perc.pca=NULL, ..., dudi=NULL){

    ## FIRST CHECKS
    grp <- as.factor(grp)
    if(length(grp) != nrow(x)) stop("Inconsistent length for grp")
    pca.select <- match.arg(pca.select)
    if(!is.null(perc.pca) & is.null(n.pca)) pca.select <- "percVar"
    if(is.null(perc.pca) & !is.null(n.pca)) pca.select <- "nbEig"
    if(!is.null(dudi) && !inherits(dudi, "dudi")) stop("dudi provided, but not of class 'dudi'")


    ## SOME GENERAL VARIABLES
    N <- nrow(x)
    REDUCEDIM <- is.null(dudi)

    if(REDUCEDIM){ # if no dudi provided
        ## PERFORM PCA ##
        maxRank <- min(dim(x))
        pcaX <- dudi.pca(x, center = center, scale = scale, scannf = FALSE, nf=maxRank)
    } else { # else use the provided dudi
        pcaX <- dudi
    }
    cumVar <- 100 * cumsum(pcaX$eig)/sum(pcaX$eig)

    if(!REDUCEDIM){
        myCol <- rep(c("black", "lightgrey"), c(ncol(pcaX$li),length(pcaX$eig)))
    } else {
        myCol <- "black"
    }

    ## select the number of retained PC for PCA
    if(is.null(n.pca) & pca.select=="nbEig"){
        plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA", col=myCol)
        cat("Choose the number PCs to retain (>=1): ")
        n.pca <- as.integer(readLines(n = 1))
    }

    if(is.null(perc.pca) & pca.select=="percVar"){
        plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA", col=myCol)
        cat("Choose the percentage of variance to retain (0-100): ")
        nperc.pca <- as.numeric(readLines(n = 1))
    }

    ## get n.pca from the % of variance to conserve
    if(!is.null(perc.pca)){
        n.pca <- min(which(cumVar >= perc.pca))
        if(perc.pca > 99.999) n.pca <- length(pcaX$eig)
        if(n.pca<1) n.pca <- 1
    }


    ## keep relevant PCs - stored in XU
    X.rank <- sum(pcaX$eig > 1e-14)
    n.pca <- min(X.rank, n.pca)
    if(n.pca >= N) n.pca <- N-1
    if(n.pca > N/3) warning("number of retained PCs of PCA may be too large (> N /3)\n results may be unstable ")
    n.pca <- round(n.pca)

    U <- pcaX$c1[, 1:n.pca, drop=FALSE] # principal axes
    rownames(U) <- colnames(x) # force to restore names
    XU <- pcaX$li[, 1:n.pca, drop=FALSE] # principal components
    XU.lambda <- sum(pcaX$eig[1:n.pca])/sum(pcaX$eig) # sum of retained eigenvalues
    names(U) <- paste("PCA-pa", 1:ncol(U), sep=".")
    names(XU) <- paste("PCA-pc", 1:ncol(XU), sep=".")


    ## PERFORM DA ##
    ldaX <- lda(XU, grp, tol=1e-30) # tol=1e-30 is a kludge, but a safe (?) one to avoid fancy rescaling by lda.default
    lda.dim <- sum(ldaX$svd^2 > 1e-10)
    ldaX$svd <- ldaX$svd[1:lda.dim]
    ldaX$scaling <- ldaX$scaling[,1:lda.dim,drop=FALSE]

    if(is.null(n.da)){
        barplot(ldaX$svd^2, xlab="Linear Discriminants", ylab="F-statistic", main="Discriminant analysis eigenvalues", col=heat.colors(length(levels(grp))) )
        cat("Choose the number discriminant functions to retain (>=1): ")
        n.da <- as.integer(readLines(n = 1))
    }

    ##n.da <- min(n.da, length(levels(grp))-1, n.pca) # can't be more than K-1 disc. func., or more than n.pca
    n.da <- round(min(n.da, lda.dim)) # can't be more than K-1 disc. func., or more than n.pca
    predX <- predict(ldaX, dimen=n.da)


    ## BUILD RESULT
    res <- list()
    res$n.pca <- n.pca
    res$n.da <- n.da
    res$tab <- XU
    res$grp <- grp
    res$var <- XU.lambda
    res$eig <- ldaX$svd^2
    res$loadings <- ldaX$scaling[, 1:n.da, drop=FALSE]
    res$means <- ldaX$means
    res$ind.coord <-predX$x
    res$grp.coord <- apply(res$ind.coord, 2, tapply, grp, mean)
    res$prior <- ldaX$prior
    res$posterior <- predX$posterior
    res$assign <- predX$class
    res$call <- match.call()


    ## optional: store loadings of variables
    if(pca.info){
        res$pca.loadings <- as.matrix(U)
        res$pca.cent <- pcaX$cent
        res$pca.norm <- pcaX$norm
        res$pca.eig <- pcaX$eig
    }

    ## optional: get loadings of variables
    if(var.contrib){
        res$var.contr <- as.matrix(U) %*% as.matrix(ldaX$scaling[,1:n.da,drop=FALSE])
        f1 <- function(x){
            temp <- sum(x*x)
            if(temp < 1e-12) return(rep(0, length(x)))
            return(x*x / temp)
        }
        res$var.contr <- apply(res$var.contr, 2, f1)
    }

    class(res) <- "dapc"
    return(res)
} # end dapc.data.frame





#############
## dapc.matrix
#############
dapc.matrix <- function(x, ...){
    return(dapc(as.data.frame(x), ...))
}




#############
## dapc.genind
#############
dapc.genind <- function(x, pop=NULL, n.pca=NULL, n.da=NULL,
                        scale=FALSE, truenames=TRUE, var.contrib=TRUE, pca.info=TRUE,
                        pca.select=c("nbEig","percVar"), perc.pca=NULL, ...){

    ## FIRST CHECKS
    if(!is.genind(x)) stop("x must be a genind object.")

    if(is.null(pop)) {
        pop.fac <- pop(x)
    } else {
        pop.fac <- pop
    }

    if(is.null(pop.fac)) stop("x does not include pre-defined populations, and `pop' is not provided")


    ## SOME GENERAL VARIABLES
    N <- nrow(x@tab)

    ## PERFORM PCA ##
    maxRank <- min(dim(x@tab))

    X <- scaleGen(x, center = TRUE, scale = scale,
                  missing = "mean", truenames = truenames)

    ## CALL DATA.FRAME METHOD ##
    res <- dapc(X, grp=pop.fac, n.pca=n.pca, n.da=n.da,
                center=FALSE, scale=FALSE, var.contrib=var.contrib,
                pca.select=pca.select, perc.pca=perc.pca)

    res$call <- match.call()

    ## restore centring/scaling
    res$pca.cent <- attr(X, "scaled:center")

    if(scale) {
        res$pca.norm <- attr(X, "scaled:scale")
    }

    return(res)
} # end dapc.genind






######################
## Function dapc.dudi
######################
dapc.dudi <- function(x, grp, ...){
    return(dapc.data.frame(x$li, grp, dudi=x, ...))
}





#################
## dapc.genlight
#################
dapc.genlight <- function(x, pop=NULL, n.pca=NULL, n.da=NULL,
                          scale=FALSE,  var.contrib=TRUE, pca.info=TRUE,
                          pca.select=c("nbEig","percVar"), perc.pca=NULL, glPca=NULL, ...){
    ## FIRST CHECKS ##
    if(!inherits(x, "genlight")) stop("x must be a genlight object.")

    pca.select <- match.arg(pca.select)

    if(is.null(pop)) {
        pop.fac <- pop(x)
    } else {
        pop.fac <- pop
    }

    if(is.null(pop.fac)) stop("x does not include pre-defined populations, and `pop' is not provided")



    ## PERFORM PCA ##
    REDUCEDIM <- is.null(glPca)

    if(REDUCEDIM){ # if no glPca provided
        maxRank <- min(c(nInd(x), nLoc(x)))
        pcaX <- glPca(x, center = TRUE, scale = scale, nf=maxRank, loadings=FALSE, returnDotProd = TRUE, ...)
    }

    if(!REDUCEDIM){ # else use the provided glPca object
        if(is.null(glPca$loadings) & var.contrib) {
            warning("Contribution of variables requested but glPca object provided without loadings.")
            var.contrib <- FALSE
        }
        pcaX <- glPca
    }

    if(is.null(n.pca)){
        cumVar <- 100 * cumsum(pcaX$eig)/sum(pcaX$eig)
    }


    ## select the number of retained PC for PCA
    if(!REDUCEDIM){
        myCol <- rep(c("black", "lightgrey"), c(ncol(pcaX$scores),length(pcaX$eig)))
    } else {
        myCol <- "black"
    }

    if(is.null(n.pca) & pca.select=="nbEig"){
        plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA", col=myCol)
        cat("Choose the number PCs to retain (>=1): ")
        n.pca <- as.integer(readLines(n = 1))
    }

    if(is.null(perc.pca) & pca.select=="percVar"){
        plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA", col=myCol)
        cat("Choose the percentage of variance to retain (0-100): ")
        nperc.pca <- as.numeric(readLines(n = 1))
    }

    ## get n.pca from the % of variance to conserve
    if(!is.null(perc.pca)){
        n.pca <- min(which(cumVar >= perc.pca))
        if(perc.pca > 99.999) n.pca <- length(pcaX$eig)
        if(n.pca<1) n.pca <- 1
    }

    if(!REDUCEDIM){
        if(n.pca > ncol(pcaX$scores)) {
            n.pca <- ncol(pcaX$scores)
        }
    }


    ## recompute PCA with loadings if needed
    if(REDUCEDIM){
        pcaX <- glPca(x, center = TRUE, scale = scale, nf=n.pca, loadings=var.contrib, matDotProd = pcaX$dotProd)
    }


    ## keep relevant PCs - stored in XU
    N <- nInd(x)
    X.rank <- sum(pcaX$eig > 1e-14)
    n.pca <- min(X.rank, n.pca)
    if(n.pca >= N) n.pca <- N-1
    if(n.pca > N/3) warning("number of retained PCs of PCA may be too large (> N /3)\n results may be unstable ")

    U <- pcaX$loadings[, 1:n.pca, drop=FALSE] # principal axes
    XU <- pcaX$scores[, 1:n.pca, drop=FALSE] # principal components
    XU.lambda <- sum(pcaX$eig[1:n.pca])/sum(pcaX$eig) # sum of retained eigenvalues
    names(U) <- paste("PCA-pa", 1:ncol(U), sep=".")
    names(XU) <- paste("PCA-pc", 1:ncol(XU), sep=".")


    ## PERFORM DA ##
    ldaX <- lda(XU, pop.fac, tol=1e-30) # tol=1e-30 is a kludge, but a safe (?) one to avoid fancy rescaling by lda.default
    lda.dim <- sum(ldaX$svd^2 > 1e-10)
    ldaX$svd <- ldaX$svd[1:lda.dim]
    ldaX$scaling <- ldaX$scaling[,1:lda.dim,drop=FALSE]

    if(is.null(n.da)){
        barplot(ldaX$svd^2, xlab="Linear Discriminants", ylab="F-statistic", main="Discriminant analysis eigenvalues", col=heat.colors(length(levels(pop.fac))) )
        cat("Choose the number discriminant functions to retain (>=1): ")
        n.da <- as.integer(readLines(n = 1))
    }

    n.da <- min(n.da, length(levels(pop.fac))-1, n.pca, sum(ldaX$svd>1e-10)) # can't be more than K-1 disc. func., or more than n.pca
    n.da <- round(n.da)
    predX <- predict(ldaX, dimen=n.da)


    ## BUILD RESULT
    res <- list()
    res$n.pca <- n.pca
    res$n.da <- n.da
    res$tab <- XU
    res$grp <- pop.fac
    res$var <- XU.lambda
    res$eig <- ldaX$svd^2
    res$loadings <- ldaX$scaling[, 1:n.da, drop=FALSE]
    res$means <- ldaX$means
    res$ind.coord <-predX$x
    res$grp.coord <- apply(res$ind.coord, 2, tapply, pop.fac, mean)
    res$prior <- ldaX$prior
    res$posterior <- predX$posterior
    res$assign <- predX$class
    res$call <- match.call()


    ## optional: store loadings of variables
    if(pca.info){
        res$pca.loadings <- as.matrix(U)
        res$pca.cent <- glMean(x,alleleAsUnit=FALSE)
        if(scale) {
            res$pca.norm <- sqrt(glVar(x,alleleAsUnit=FALSE))
        } else {
            res$pca.norm <- rep(1, nLoc(x))
        }
        res$pca.eig <- pcaX$eig
    }

    ## optional: get loadings of variables
    if(var.contrib){
        res$var.contr <- as.matrix(U) %*% as.matrix(ldaX$scaling[,1:n.da,drop=FALSE])
        f1 <- function(x){
            temp <- sum(x*x)
            if(temp < 1e-12) return(rep(0, length(x)))
            return(x*x / temp)
        }
        res$var.contr <- apply(res$var.contr, 2, f1)
    }

    class(res) <- "dapc"
    return(res)
} # end dapc.genlight






######################
# Function print.dapc
######################
print.dapc <- function(x, ...){
    cat("\t#################################################\n")
    cat("\t# Discriminant Analysis of Principal Components #\n")
    cat("\t#################################################\n")
    cat("class: ")
    cat(class(x))
    cat("\n$call: ")
    print(x$call)
    cat("\n$n.pca:", x$n.pca, "first PCs of PCA used")
    cat("\n$n.da:", x$n.da, "discriminant functions saved")
    cat("\n$var (proportion of conserved variance):", round(x$var,3))
    cat("\n\n$eig (eigenvalues): ")
    l0 <- sum(x$eig >= 0)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5)
        cat(" ...\n\n")

    ## vectors
    TABDIM <- 4
    if(!is.null(x$pca.loadings)){
        TABDIM <- TABDIM + 3
    }
    sumry <- array("", c(TABDIM, 3), list(1:TABDIM, c("vector", "length", "content")))
    sumry[1, ] <- c('$eig', length(x$eig),  'eigenvalues')
    sumry[2, ] <- c('$grp', length(x$grp), 'prior group assignment')
    sumry[3, ] <- c('$prior', length(x$prior), 'prior group probabilities')
    sumry[4, ] <- c('$assign', length(x$assign), 'posterior group assignment')
    if(!is.null(x$pca.loadings)){
        sumry[5, ] <- c('$pca.cent', length(x$pca.cent), 'centring vector of PCA')
        sumry[6, ] <- c('$pca.norm', length(x$pca.norm), 'scaling vector of PCA')
        sumry[7, ] <- c('$pca.eig', length(x$pca.eig), 'eigenvalues of PCA')
    }
    class(sumry) <- "table"
    print(sumry)

    ## data.frames
    cat("\n")
    TABDIM <- 6
    if(!is.null(x$pca.loadings)){
        TABDIM <- TABDIM + 1
    }
    if(!is.null(x$var.contr)){
        TABDIM <- TABDIM + 1
    }

    sumry <- array("", c(TABDIM, 4), list(1:TABDIM, c("data.frame", "nrow", "ncol", "content")))

    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "retained PCs of PCA")
    sumry[2, ] <- c("$means", nrow(x$means), ncol(x$means), "group means")
    sumry[3, ] <- c("$loadings", nrow(x$loadings), ncol(x$loadings), "loadings of variables")
    sumry[4, ] <- c("$ind.coord", nrow(x$ind.coord), ncol(x$ind.coord), "coordinates of individuals (principal components)")
    sumry[5, ] <- c("$grp.coord", nrow(x$grp.coord), ncol(x$grp.coord), "coordinates of groups")
    sumry[6, ] <- c("$posterior", nrow(x$posterior), ncol(x$posterior), "posterior membership probabilities")
    if(!is.null(x$pca.loadings)){
        sumry[7, ] <- c("$pca.loadings", nrow(x$pca.loadings), ncol(x$pca.loadings), "PCA loadings of original variables")
    }
    if(!is.null(x$var.contr)){
        sumry[TABDIM, ] <- c("$var.contr", nrow(x$var.contr), ncol(x$var.contr), "contribution of original variables")
    }
    class(sumry) <- "table"
    print(sumry)

    ## cat("\nother elements: ")
    ## if (length(names(x)) > 15)
    ##     cat(names(x)[15:(length(names(x)))], "\n")
    ## else cat("NULL\n")
    cat("\n")
} # end print.dapc






##############
## summary.dapc
##############
summary.dapc <- function(object, ...){
    x <- object
    res <- list()

    ## number of dimensions
    res$n.dim <- ncol(x$loadings)
    res$n.pop <- length(levels(x$grp))

    ## assignment success
    temp <- as.character(x$grp)==as.character(x$assign)
    res$assign.prop <- mean(temp)
    res$assign.per.pop <- tapply(temp, x$grp, mean)

    ## group sizes
    res$prior.grp.size <- table(x$grp)
    res$post.grp.size <- table(x$assign)

    return(res)
} # end summary.dapc






##############
## scatter.dapc
##############
scatter.dapc <- function(x, xax=1, yax=2, grp=x$grp, col=seasun(length(levels(grp))), pch=20, bg="white", solid=.7,
                         scree.da=TRUE, scree.pca=FALSE, posi.da="bottomright", posi.pca="bottomleft", bg.inset="white",
                         ratio.da=.25, ratio.pca=.25, inset.da=0.02, inset.pca=0.02, inset.solid=.5,
                         onedim.filled=TRUE, mstree=FALSE, lwd=1, lty=1, segcol="black",
                         legend=FALSE, posi.leg="topright", cleg=1, txt.leg=levels(grp),
                         cstar = 1, cellipse = 1.5, axesell = FALSE, label = levels(grp), clabel = 1, xlim = NULL, ylim = NULL,
                         grid = FALSE, addaxes = TRUE, origin = c(0,0), include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft",
                         cgrid = 1, pixmap = NULL, contour = NULL, area = NULL, ...){
    ONEDIM <- xax==yax | ncol(x$ind.coord)==1


    ## recycle color and pch
    col <- rep(col, length(levels(grp)))
    pch <- rep(pch, length(levels(grp)))
    col <- transp(col, solid)
    bg.inset <- transp(bg.inset, inset.solid)

    ## handle grp
    if(is.null(grp)){
        grp <- x$grp
    }

    ## handle xax or yax NULL
    if(is.null(xax)||is.null(yax)){
        xax <- 1L
        yax <- ifelse(ncol(x$ind.coord)==1L, 1L, 2L)
        ONEDIM <- TRUE
    }

    ## handle 1 dimensional plot
    if(!ONEDIM){
        ## set par
        opar <- par(mar = par("mar"))
        par(mar = c(0.1, 0.1, 0.1, 0.1), bg=bg)
        on.exit(par(opar))
        axes <- c(xax,yax)

        ## basic empty plot
        s.class(x$ind.coord[,axes], fac=grp, col=col, cpoint=0, cstar = cstar, cellipse = cellipse, axesell = axesell, label = label,
                clabel = clabel, xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, origin = origin, include.origin = include.origin,
                sub = sub, csub = csub, possub = possub, cgrid = cgrid, pixmap = pixmap, contour = contour, area = area)

        ## add points
        colfac <- pchfac <- grp
        levels(colfac) <- col
        levels(pchfac) <- pch
        colfac <- as.character(colfac)
        pchfac <- as.character(pchfac)
        if(is.numeric(col)) colfac <- as.numeric(colfac)
        if(is.numeric(pch)) pchfac <- as.numeric(pchfac)

        points(x$ind.coord[,xax], x$ind.coord[,yax], col=colfac, pch=pchfac, ...)
        s.class(x$ind.coord[,axes], fac=grp, col=col, cpoint=0, add.plot=TRUE, cstar = cstar, cellipse = cellipse, axesell = axesell, label = label,
                clabel = clabel, xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, origin = origin, include.origin = include.origin,
                sub = sub, csub = csub, possub = possub, cgrid = cgrid, pixmap = pixmap, contour = contour, area = area)

        ## add minimum spanning tree if needed
        if(mstree){
            meanposi <- apply(x$tab,2, tapply, grp, mean)
            D <- dist(meanposi)^2
            tre <- ade4::mstree(D)
            x0 <- x$grp.coord[tre[,1], axes[1]]
            y0 <- x$grp.coord[tre[,1], axes[2]]
            x1 <- x$grp.coord[tre[,2], axes[1]]
            y1 <- x$grp.coord[tre[,2], axes[2]]
            segments(x0, y0, x1, y1, lwd=lwd, lty=lty, col=segcol)
        }

    } else {
        ## set screeplot of DA to FALSE (just 1 bar)
        scree.da <- FALSE

        ## get plotted axis
        if(ncol(x$ind.coord)==1) {
            pcLab <- 1
        } else{
            pcLab <- xax
        }
        ## get densities
        ldens <- tapply(x$ind.coord[,pcLab], grp, density)
        allx <- unlist(lapply(ldens, function(e) e$x))
        ally <- unlist(lapply(ldens, function(e) e$y))
        par(bg=bg)
        plot(allx, ally, type="n", xlab=paste("Discriminant function", pcLab), ylab="Density")
        for(i in 1:length(ldens)){
            if(!onedim.filled) {
                lines(ldens[[i]]$x,ldens[[i]]$y, col=col[i], lwd=2) # add lines
            } else {
                polygon(c(ldens[[i]]$x,rev(ldens[[i]]$x)),c(ldens[[i]]$y,rep(0,length(ldens[[i]]$x))), col=col[i], lwd=2, border=col[i]) # add lines
            }
            points(x=x$ind.coord[grp==levels(grp)[i],pcLab], y=rep(0, sum(grp==levels(grp)[i])), pch="|", col=col[i]) # add points for indiv
        }
    }

    ## ADD INSETS ##
    ## group legend
    if(legend){
        ## add a legend
        temp <- list(...)$cex
        if(is.null(temp)) temp <- 1
        if(ONEDIM | temp<0.5 | all(pch=="")) {
            legend(posi.leg, fill=col, legend=txt.leg, cex=cleg, bg=bg.inset)
        } else {
            legend(posi.leg, col=col, legend=txt.leg, cex=cleg, bg=bg.inset, pch=pch, pt.cex=temp)
        }
    }

    ## eigenvalues discriminant analysis
    if(scree.da && ratio.da>.01) {
        inset <- function(){
            myCol <- rep("white", length(x$eig))
            myCol[1:x$n.da] <- "grey"
            myCol[c(xax, yax)] <- "black"
            myCol <- transp(myCol, inset.solid)
            barplot(x$eig, col=myCol, xaxt="n", yaxt="n", ylim=c(0, x$eig[1]*1.1))
            mtext(side=3, "DA eigenvalues", line=-1.2, adj=.8)
            box()
        }

        add.scatter(inset(), posi=posi.da, ratio=ratio.da, bg.col=bg.inset, inset=inset.da)
        ##add.scatter.eig(x$eig, ncol(x$loadings), axes[1], axes[2], posi=posi, ratio=ratio, csub=csub) # does not allow for bg
    }

    ## eigenvalues PCA
    if(scree.pca && !is.null(x$pca.eig) && ratio.pca>.01) {
        inset <- function(){
            temp <- 100* cumsum(x$pca.eig) / sum(x$pca.eig)
            myCol <- rep(c("black","grey"), c(x$n.pca, length(x$pca.eig)))
            myCol <- transp(myCol, inset.solid)
            plot(temp, col=myCol, ylim=c(0,115),
                 type="h", xaxt="n", yaxt="n", xlab="", ylab="", lwd=2)
            mtext(side=3, "PCA eigenvalues", line=-1.2, adj=.1)
        }
        add.scatter(inset(), posi=posi.pca, ratio=ratio.pca, bg.col=bg.inset, inset=inset.pca)
    }


    return(invisible(match.call()))
} # end scatter.dapc






############
## assignplot
############
assignplot <- function(x, only.grp=NULL, subset=NULL, new.pred=NULL, cex.lab=.75, pch=3){
    if(!inherits(x, "dapc")) stop("x is not a dapc object")

    ## handle data from predict.dapc ##
    if(!is.null(new.pred)){
        n.new <- length(new.pred$assign)
        x$grp <- c(as.character(x$grp), rep("unknown", n.new))
        x$assign <- c(as.character(x$assign), as.character(new.pred$assign))
        x$posterior <- rbind(x$posterior, new.pred$posterior)
    }


    ## treat other arguments ##
    if(!is.null(only.grp)){
        only.grp <- as.character(only.grp)
        ori.grp <- as.character(x$grp)
        x$grp <- x$grp[only.grp==ori.grp]
        x$assign <- x$assign[only.grp==ori.grp]
        x$posterior <- x$posterior[only.grp==ori.grp, , drop=FALSE]
    } else if(!is.null(subset)){
        x$grp <- x$grp[subset]
        x$assign <- x$assign[subset]
        x$posterior <- x$posterior[subset, , drop=FALSE]
    }


    ##table.paint(x$posterior, col.lab=ori.grp, ...)
    ## symbols(x$posterior)


    ## FIND PLOT PARAMETERS
    n.grp <- ncol(x$posterior)
    n.ind <- nrow(x$posterior)
    Z <- t(x$posterior)
    Z <- Z[,ncol(Z):1,drop=FALSE ]

    image(x=1:n.grp, y=seq(.5, by=1, le=n.ind), Z, col=rev(heat.colors(100)), yaxt="n", ylab="", xaxt="n", xlab="Clusters")
    axis(side=1, at=1:n.grp,tick=FALSE, labels=colnames(x$posterior))
    axis(side=2, at=seq(.5, by=1, le=n.ind), labels=rev(rownames(x$posterior)), las=1, cex.axis=cex.lab)
    abline(h=1:n.ind, col="lightgrey")
    abline(v=seq(0.5, by=1, le=n.grp))
    box()

    newGrp <- colnames(x$posterior)
    x.real.coord <- rev(match(x$grp, newGrp))
    y.real.coord <- seq(.5, by=1, le=n.ind)

    points(x.real.coord, y.real.coord, col="deepskyblue2", pch=pch)

    return(invisible(match.call()))
} # end assignplot





############
## compoplot
############
compoplot <- function(x, only.grp=NULL, subset=NULL, new.pred=NULL, col=NULL, lab=NULL,
                      legend=TRUE, txt.leg=NULL, ncol=4, posi=NULL, cleg=.8, bg=transp("white"), ...){
    if(!inherits(x, "dapc")) stop("x is not a dapc object")


    ## HANDLE ARGUMENTS ##
    ngrp <- length(levels(x$grp))

    ## col
    if(is.null(col)){
        col <- rainbow(ngrp)
    }

    ## lab
    if(is.null(lab)){
        lab <- rownames(x$tab)
    } else {
        ## recycle labels
       lab <- rep(lab, le=nrow(x$tab))
    }

    ## posi
    if(is.null(posi)){
        posi <- list(x=0, y=-.01)
    }

    ## txt.leg
    if(is.null(txt.leg)){
        txt.leg <- levels(x$grp)
    }

    ## HANDLE DATA FROM PREDICT.DAPC ##
    if(!is.null(new.pred)){
        n.new <- length(new.pred$assign)
        x$grp <- c(as.character(x$grp), rep("unknown", n.new))
        x$assign <- c(as.character(x$assign), as.character(new.pred$assign))
        x$posterior <- rbind(x$posterior, new.pred$posterior)
        lab <- c(lab, rownames(new.pred$posterior))
    }


    ## TREAT OTHER ARGUMENTS ##
    if(!is.null(only.grp)){
        only.grp <- as.character(only.grp)
        ori.grp <- as.character(x$grp)
        x$grp <- x$grp[only.grp==ori.grp]
        x$assign <- x$assign[only.grp==ori.grp]
        x$posterior <- x$posterior[only.grp==ori.grp, , drop=FALSE]
        lab <- lab[only.grp==ori.grp]
    } else if(!is.null(subset)){
        x$grp <- x$grp[subset]
        x$assign <- x$assign[subset]
        x$posterior <- x$posterior[subset, , drop=FALSE]
        lab <- lab[subset]
    }


    ## MAKE THE PLOT ##
    Z <- t(x$posterior)
    barplot(Z, border=NA, col=col, ylab="membership probability", names=lab, las=3, ...)

    if(legend){
        oxpd <- par("xpd")
        par(xpd=TRUE)
        legend(posi, fill=col, leg=txt.leg, cex=cleg, ncol=ncol, bg=bg)
        on.exit(par(xpd=oxpd))
    }

    return(invisible(match.call()))
} # end compoplot





###############
## a.score
###############
a.score <- function(x, n.sim=10, ...){
    if(!inherits(x,"dapc")) stop("x is not a dapc object")

    ## perform DAPC based on permuted groups
    lsim <- lapply(1:n.sim, function(i) summary(dapc(x$tab, sample(x$grp), n.pca=x$n.pca, n.da=x$n.da))$assign.per.pop)
    sumry <- summary(x)

    ## get the a-scores
    f1 <- function(Pt, Pf){
        tol <- 1e-7
        ##res <- (Pt-Pf) / (1-Pf)
        ##res[Pf > (1-tol)] <- 0
        res <- Pt-Pf
        return(res)
    }

    lscores <- lapply(lsim, function(e) f1(sumry$assign.per.pop, e))

    ## make a table of a-scores
    tab <- data.frame(lscores)
    colnames(tab) <- paste("sim", 1:n.sim, sep=".")
    rownames(tab) <- names(sumry$assign.per.pop)
    tab <- t(as.matrix(tab))

    ## make result
    res <- list()
    res$tab <- tab
    res$pop.score <- apply(tab, 2, mean)
    res$mean <- mean(tab)

    return(res)

} # end a.score







##############
## optim.a.score
##############
optim.a.score <- function(x, n.pca=1:ncol(x$tab), smart=TRUE, n=10, plot=TRUE,
                         n.sim=10, n.da=length(levels(x$grp)), ...){
    ## A FEW CHECKS ##
    if(!inherits(x,"dapc")) stop("x is not a dapc object")
    if(max(n.pca)>ncol(x$tab)) {
        n.pca <- min(n.pca):ncol(x$tab)
    }
    if(n.da>length(levels(x$grp))){
        n.da <- min(n.da):length(levels(x$grp))
    }
    pred <- NULL
    if(length(n.pca)==1){
        n.pca <- 1:n.pca
    }
    if(length(n.da)==1){
        n.da <- 1:n.da
    }


    ## AUXILIARY FUNCTION ##
    f1 <- function(ndim){
        temp <- dapc(x$tab[,1:ndim,drop=FALSE], x$grp, n.pca=ndim, n.da=x$n.da)
        a.score(temp, n.sim=n.sim)$pop.score
    }


    ## SMART: COMPUTE A FEW VALUES, PREDICT THE BEST PICK ##
    if(smart){
        if(!require(stats)) stop("the package stats is required for 'smart' option")
        o.min <- min(n.pca)
        o.max <- max(n.pca)
        n.pca <- pretty(n.pca, n) # get evenly spaced nb of retained PCs
        n.pca <- n.pca[n.pca>0 & n.pca<=ncol(x$tab)]
        if(!any(o.min==n.pca)) n.pca <- c(o.min, n.pca) # make sure range is OK
        if(!any(o.max==n.pca)) n.pca <- c(o.max, n.pca) # make sure range is OK
        lres <- lapply(n.pca, f1)
        names(lres) <- n.pca
        means <- sapply(lres, mean)
        sp1 <- smooth.spline(n.pca, means) # spline smoothing
        pred <- predict(sp1, x=1:max(n.pca))
        best <- pred$x[which.max(pred$y)]
    } else { ## DO NOT TRY TO BE SMART ##
        lres <- lapply(n.pca, f1)
        names(lres) <- n.pca
        best <- which.max(sapply(lres, mean))
        means <- sapply(lres, mean)
    }


    ## MAKE FINAL OUTPUT ##
    res <- list()
    res$pop.score <- lres
    res$mean <- means
    if(!is.null(pred)) res$pred <- pred
    res$best <- best

    ## PLOTTING (OPTIONAL) ##
    if(plot){
        if(smart){
            boxplot(lres, at=n.pca, col="gold", xlab="Number of retained PCs", ylab="a-score", xlim=range(n.pca)+c(-1,1), ylim=c(-.1,1.1))
            lines(pred, lwd=3)
            points(pred$x[best], pred$y[best], col="red", lwd=3)
            title("a-score optimisation - spline interpolation")
            mtext(paste("Optimal number of PCs:", res$best), side=3)
        } else {
            myCol <- rep("gold", length(lres))
            myCol[best] <- "red"
            boxplot(lres, at=n.pca, col=myCol, xlab="Number of retained PCs", ylab="a-score", xlim=range(n.pca)+c(-1,1), ylim=c(-.1,1.1))
            lines(n.pca, sapply(lres, mean), lwd=3, type="b")
            myCol <- rep("black", length(lres))
            myCol[best] <- "red"
            points(n.pca, res$mean, lwd=3, col=myCol)
            title("a-score optimisation - basic search")
            mtext(paste("Optimal number of PCs:", res$best), side=3)
        }
    }

    return(res)
} # end optim.a.score






#############
## as.lda.dapc
#############
as.lda <- function(...){
    UseMethod("as.lda")
}

as.lda.dapc <- function(x, ...){
    if(!inherits(x,"dapc")) stop("x is not a dapc object")
    res <- list()

    res$N <- nrow(res$ind.coord)
    res$call <- match.call()
    res$counts <- as.integer(table(x$grp))
    res$lev <- names(res$counts) <- levels(x$grp)
    res$means <- x$means
    res$prior <- x$prior
    res$scaling <- x$loadings
    res$svd <- sqrt(x$eig)

    class(res) <- "lda"

    return(res)
} # end as.lda.dapc






##############
## predict.dapc
##############
predict.dapc <- function(object, newdata, prior = object$prior, dimen,
                         method = c("plug-in", "predictive", "debiased"), ...){

    if(!inherits(object,"dapc")) stop("x is not a dapc object")
    method <- match.arg(method)

    x <- as.lda(object)


    ## HANDLE NEW DATA ##
    if(!missing(newdata)){
        ## make a few checks
        if(is.null(object$pca.loadings)) stop("DAPC object does not contain loadings of original variables. \nPlease re-run DAPC using 'pca.loadings=TRUE'.")
        newdata <- as.matrix(newdata) # to force conversion, notably from genlight objects
        if(ncol(newdata) != nrow(object$pca.loadings)) stop("Number of variables in newdata does not match original data.")

        ## centre/scale data
        for(i in 1:nrow(newdata)){ # this is faster for large, flat matrices)
            newdata[i,] <- (newdata[i,] - object$pca.cent) / object$pca.norm
        }
        newdata[is.na(newdata)] <- 0

        ## project as supplementary individuals
        XU <- newdata %*% as.matrix(object$pca.loadings)
    } else {
        XU <- object$tab
    }

    ## FORCE IDENTICAL VARIABLE NAMES ##
    colnames(XU) <- colnames(object$tab)


    ## HANDLE DIMEN ##
    if(!missing(dimen)){
        if(dimen > object$n.da) stop(paste("Too many dimensions requested. \nOnly", object$n.da, "discriminant functions were saved in DAPC."))
    } else {
        dimen <- object$n.da
    }

    ## CALL PREDICT.LDA ##
    temp <- predict(x, XU, prior, dimen, method, ...)


    ## FORMAT OUTPUT ##
    res <- list()
    res$assign <- temp$class
    res$posterior <- temp$posterior
    res$ind.scores <- temp$x

    return(res)

} # end predict.dapc







############
## xvalDapc
############

xvalDapc <- function (x, ...) UseMethod("xvalDapc")

xvalDapc.data.frame <- function(x, grp, n.pca.max, n.da=NULL, training.set = 0.9,
                                result=c("groupMean","overall"),
                                center=TRUE, scale=FALSE, n.pca=NULL, n.rep=10, ...){

    ## CHECKS ##
    grp <- factor(grp)
    n.pca <- n.pca[n.pca>0]
    result <- match.arg(result)
    if(is.null(n.da)) {
        n.da <- length(levels(grp))-1
    }

    ## GET TRAINING SET SIZE ##
    N <- nrow(x)
    N.training <- round(N*training.set)

    ## GET FULL PCA ##
    if(missing(n.pca.max)) n.pca.max <- min(dim(x))
    pcaX <- dudi.pca(x, nf=n.pca.max, scannf=FALSE, center=center, scale=scale)
    n.pca.max <- min(n.pca.max,pcaX$rank,N.training-1)

    ## DETERMINE N.PCA IF NEEDED ##
    if(is.null(n.pca)){
        n.pca <- round(pretty(1:n.pca.max,10))
    }
    n.pca <- n.pca[n.pca>0 & n.pca<(N.training-1)]

    ## FUNCTION GETTING THE % OF ACCURATE PREDICTION FOR ONE NUMBER OF PCA PCs ##
    ## n.pca is a number of retained PCA PCs
    VOID.GRP <- FALSE # will be TRUE if empty group happened
    get.prop.pred <- function(n.pca){
        f1 <- function(){
            toKeep <- sample(1:N, N.training)
            if(!(all(table(grp[toKeep])>0) & all(table(grp[-toKeep])>0))) VOID.GRP <<- TRUE
            temp.pca <- pcaX
            temp.pca$li <- temp.pca$li[toKeep,,drop=FALSE]
            temp.dapc <- suppressWarnings(dapc(x[toKeep,,drop=FALSE], grp[toKeep], n.pca=n.pca, n.da=n.da, dudi=temp.pca))
            temp.pred <- predict.dapc(temp.dapc, newdata=x[-toKeep,,drop=FALSE])
            if(result=="overall"){
                out <- mean(temp.pred$assign==grp[-toKeep])
            }
            if(result=="groupMean"){
                out <- mean(tapply(temp.pred$assign==grp[-toKeep], grp[-toKeep], mean), na.rm=TRUE)
            }
            return(out)
        }
        return(replicate(n.rep, f1()))
    }


    ## GET %SUCCESSFUL OF ACCURATE PREDICTION FOR ALL VALUES ##
    res.all <- unlist(lapply(n.pca, get.prop.pred))
    if(VOID.GRP) warning("At least one group was absent from the training / validating sets.\nTry using smaller training sets.")
    res <- data.frame(n.pca=rep(n.pca, each=n.rep), success=res.all)
    return(res)
} # end xvalDapc.data.frame


xvalDapc.matrix <- xvalDapc.data.frame






## #############
## ## discriVal
## #############

## discriVal <- function (x, ...) UseMethod("discriVal")

## discriVal.data.frame <- function(x, grp, n.pca.max, n.da=NULL, center=TRUE, scale=FALSE, n.pca=NULL, ...){

##     ## CHECKS ##
##     grp <- factor(grp)
##     n.pca <- n.pca[n.pca>0]
##     if(is.null(n.da)) {
##         n.da <- length(levels(grp))-1
##     }

##     ## GET FULL PCA ##
##     if(missing(n.pca.max)) n.pca.max <- min(dim(x))
##     pcaX <- dudi.pca(x, nf=n.pca.max, scannf=FALSE, center=center, scale=scale)
##     n.pca.max <- min(n.pca.max,pcaX$rank)

##     ## DETERMINE N.PCA IF NEEDED ##
##     if(is.null(n.pca)){
##         n.pca <- round(pretty(1:n.pca.max,10))
##     }
##     n.pca <- n.pca[n.pca>0 & n.pca<n.pca.max]

##     ## FUNCTION GETTING THE TOTAL DISCRIMINATION (SUM OF EIGENVALUES) FOR ONE GIVEN NB OF PCA PCs ##
##     ## n.pca is a number of retained PCA PCs
##     get.totdiscr <- function(n.pca){
##             temp.dapc <- suppressWarnings(dapc(x, grp, n.pca=n.pca, n.da=n.da, dudi=pcaX))
##             return(sum(temp.dapc$eig))
##     }


##     ## GET %SUCCESSFUL OF ACCURATE PREDICTION FOR ALL VALUES ##
##     res.all <- sapply(n.pca, get.totdiscr)
##     res <- data.frame(n.pca=n.pca, success=res.all)
##     return(res)
## } # end discriVal.data.frame


## discriVal.matrix <- discriVal.data.frame




## There's a bunch of problems down there, commenting it for nowé
## xval.dapc <- function(object, n.pca, n.da, training.set = 90, ...){
##   training.set = training.set/100
##   kept.id <- unlist(tapply(1:nInd(object), pop(object), function(e) {pop.size = length(e); pop.size.train = round(pop.size * training.set); sample(e, pop.size.train, replace=FALSE)})) # this can't work: nInd/pop not defined for DAPC objects
##   training <- object[kept.id]
##   validating <- object[-kept.id]
##   post = vector(mode = 'list', length = n.pca)
##   asgn = vector(mode = 'list', length = n.pca)
##   ind = vector(mode = 'list', length = n.pca)
##   mtch = vector(mode = 'list', length = n.pca)
##   for(i in 1:n.pca){
##     dapc.base = dapc(training, n.pca = i, n.da = 15) # Why 15??
##     dapc.p = predict.dapc(dapc.base, newdata = validating)
##     match.prp = mean(as.character(dapc.p$assign)==as.character(pop(validating)))
##     post[[i]] = dapc.p$posterior
##     asgn[[i]] = dapc.p$assign
##     ind[[i]] = dapc.p$ind.score
##     mtch[[i]] = match.prp
##   }
##   res = list(assign = asgn, posterior = post, ind.score = ind, match.prp = mtch)
##   return(res)
## } # end of xval.dapc

## xval.genind  <- function(object, n.pca, n.da, training.set = 90, ...){
##   res = xval.dapc(object = object, n.pca = n.pca, n.da = n.da, training.set = training.set)
##   return(res)
## }
## ###############
## ## randtest.dapc
## ###############
## ##randtest.dapc <- function(x, nperm = 999, ...){

## ##} # end randtest.dapc




######## TESTS IN R #######

## TEST PREDICT.DAPC ##
## data(sim2pop)
## temp <- seppop(sim2pop)
## temp <- lapply(temp, function(e) hybridize(e,e,n=30)) # force equal pop sizes
## hyb <- hybridize(temp[[1]], temp[[2]], n=30)
## newdat <- repool(temp[[1]], temp[[2]], hyb)
## pop(newdat) <- rep(c("pop A", "popB", "hyb AB"), c(30,30,30))


## ##dapc1 <- dapc(newdat[1:61],n.pca=10,n.da=1)
## dapc1 <- dapc(newdat[1:60],n.pca=2,n.da=1)
## scatter(dapc1)
## hyb.pred <- predict(dapc1, newdat[61:90])

## scatter(dapc1)
## points(hyb.pred$ind.scores, rep(.1, 30))

## assignplot(dapc1, new.pred=hyb.pred)
## title("30 indiv popA, 30 indiv pop B, 30 hybrids")
