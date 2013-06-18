\encoding{UTF-8}
\name{dapc}
\alias{dapc}
\alias{dapc.data.frame}
\alias{dapc.matrix}
\alias{dapc.genind}
\alias{dapc.dudi}
\alias{dapc.genlight}
\alias{print.dapc}
\alias{summary.dapc}
\alias{predict.dapc}
\alias{as.lda}
\alias{as.lda.dapc}
\title{Discriminant Analysis of Principal Components (DAPC)}
\description{
  These functions implement the Discriminant Analysis of Principal
  Components (DAPC, Jombart et al. 2010). This method descibes the
  diversity between pre-defined groups. When groups are unknown, use
  \code{find.clusters} to infer genetic clusters. See 'details' section
  for a succint description of the method, and
  \code{vignette("adegenet-dapc")} for a tutorial. Graphical methods for
  DAPC are documented in \code{\link{scatter.dapc}} (see
  \code{?scatter.dapc}).

 \code{dapc} is a generic function performing the DAPC on the following
 types of objects:\cr
 - \code{data.frame} (only numeric data)\cr
 - \code{matrix} (only numeric data)\cr
 - \code{\linkS4class{genind}} objects (genetic markers)\cr
 - \code{\linkS4class{genlight}} objects (genome-wide SNPs)

 These methods all return an object with class \code{dapc}.

 Functions that can be applied to these objects are (the ".dapc" can be
 ommitted):

  - \code{print.dapc}: prints the content of a \code{dapc} object.\cr
  - \code{summary.dapc}: extracts useful information from a \code{dapc}
  object.\cr
  - \code{predict.dapc}: predicts group memberships based on DAPC results.\cr
  - \code{xvalDapc}: performs cross-validation of DAPC using varying
  numbers of PCs (and keeping the number of discriminant functions
  fixed); it currently has methods for \code{data.frame} and \code{matrix}.\cr



  DAPC implementation calls upon \code{\link[ade4]{dudi.pca}} from the
  \code{ade4} package (except for \linkS4class{genlight} objects)
  and \code{\link[MASS]{lda}} from the \code{MASS} package. The
  \code{predict} procedure uses \code{\link[MASS]{predict.lda}} from the
  \code{MASS} package.

  \code{as.lda} is a generic with a method for \code{dapc} object which
  converts these objects into outputs similar to that of
  \code{lda.default}.
}
\usage{
\method{dapc}{data.frame}(x, grp, n.pca=NULL, n.da=NULL, center=TRUE,
     scale=FALSE,var.contrib=TRUE, pca.info=TRUE, pca.select=c("nbEig","percVar"),
     perc.pca=NULL, \ldots, dudi=NULL)

\method{dapc}{matrix}(x, \ldots)

\method{dapc}{genind}(x, pop=NULL, n.pca=NULL, n.da=NULL, scale=FALSE,
     truenames=TRUE, var.contrib=TRUE, pca.info=TRUE,
     pca.select=c("nbEig","percVar"), perc.pca=NULL, \ldots)

\method{dapc}{genlight}(x, pop = NULL, n.pca = NULL, n.da = NULL, scale
    = FALSE, var.contrib = TRUE, pca.info=TRUE, pca.select = c("nbEig", "percVar"),
    perc.pca = NULL, glPca = NULL, \ldots)

\method{dapc}{dudi}(x, grp, \ldots)

\method{print}{dapc}(x, \dots)

\method{summary}{dapc}(object, \dots)

\method{predict}{dapc}(object, newdata, prior = object$prior, dimen,
         method = c("plug-in", "predictive", "debiased"), ...)
}
\arguments{
  \item{x}{\code{a data.frame}, \code{matrix}, or \code{\linkS4class{genind}}
    object. For the \code{data.frame} and \code{matrix} arguments, only
    quantitative variables should be provided.}
  \item{grp,pop}{a \code{factor} indicating the group membership of
    individuals; for \code{scatter}, an optional grouping of individuals.}
  \item{n.pca}{an \code{integer} indicating the number of axes retained in the
    Principal Component Analysis (PCA) step. If \code{NULL}, interactive selection is triggered.}
  \item{n.da}{an \code{integer} indicating the number of axes retained in the
    Discriminant Analysis step. If \code{NULL}, interactive selection is triggered.}
  \item{center}{a \code{logical} indicating whether variables should be centred to
    mean 0 (TRUE, default) or not (FALSE). Always TRUE for \linkS4class{genind} objects.}
  \item{scale}{a \code{logical} indicating whether variables should be scaled
    (TRUE) or not (FALSE, default). Scaling consists in dividing variables by their
    (estimated) standard deviation to account for trivial differences in
    variances.}
  \item{var.contrib}{a \code{logical} indicating whether the
    contribution of original variables (alleles, for \linkS4class{genind} objects)
    should be provided (TRUE, default) or not (FALSE). Such output can be useful,
    but can also create huge matrices when there is a lot of variables.}
  \item{pca.info}{a \code{logical} indicating whether information about
    the prior PCA should be stored (TRUE, default) or not (FALSE). This
    information is required to predict group membership of new
    individuals using \code{predict}, but makes the object slightly bigger.}
  \item{pca.select}{a \code{character} indicating the mode of selection of PCA
    axes, matching either "nbEig" or "percVar". For "nbEig", the user
    has to specify the number of axes retained (interactively, or via
    \code{n.pca}). For "percVar", the user has to specify the minimum amount of
    the total variance to be preserved by the retained axes, expressed as a
    percentage (interactively, or via \code{perc.pca}).  }
  \item{perc.pca}{a \code{numeric} value between 0 and 100 indicating the
    minimal percentage of the total variance of the data to be expressed by the
    retained axes of PCA.}
  \item{\ldots}{further arguments to be passed to other functions. For
    \code{dapc.matrix}, arguments are to match those of
    \code{dapc.data.frame}; for \code{dapc.genlight}, arguments passed
    to \code{\link{glPca}}}
  \item{glPca}{an optional \code{\link{glPca}} object; if provided,
    dimension reduction is not performed (saving computational time) but
    taken directly from this object.}
  \item{object}{a \code{dapc} object.}
  \item{truenames}{a \code{logical} indicating whether true (i.e., user-specified)
    labels should be used in object outputs (TRUE, default) or not (FALSE).}
  \item{dudi}{optionally, a multivariate analysis with the class
    \code{dudi} (from the ade4 package). If provided, prior PCA will be
    ignored, and this object will be used as a prior step for variable
    orthogonalisation.}
  \item{newdata}{an optional dataset of individuals whose membership is
    seeked; can be a data.frame, a matrix, a \linkS4class{genind} or a
    \linkS4class{genlight} object, but object class must match the
    original ('training') data. In particular, variables must be exactly
    the same as in the original data. For  \linkS4class{genind}
    objects, see \code{\link{repool}} to ensure matching of alleles.}
  \item{prior,dimen,method}{see \code{?predict.lda}.}
}
\details{
  The Discriminant Analysis of Principal Components (DAPC) is designed
  to investigate the genetic structure of biological populations. This
  multivariate method consists in a two-steps procedure. First, genetic
  data are transformed (centred, possibly scaled) and submitted to a
  Principal Component Analysis (PCA). Second, principal components of
  PCA are submitted to a Linear Discriminant Analysis (LDA). A trivial
  matrix operation allows to express discriminant functions as linear
  combination of alleles, therefore allowing one to compute allele
  contributions. More details about the computation of DAPC are to be
  found in the indicated reference.

  DAPC does not infer genetic clusters ex nihilo; for this, see the
  \code{\link{find.clusters}} function.
}
\value{
  === dapc objects ===\cr
  The class \code{dapc} is a list with the following
  components:\cr
  \item{call}{the matched call.}
  \item{n.pca}{number of PCA axes retained}
  \item{n.da}{number of DA axes retained}
  \item{var}{proportion of variance conserved by PCA principal components}
  \item{eig}{a numeric vector of eigenvalues.}
  \item{grp}{a factor giving prior group assignment}
  \item{prior}{a numeric vector giving prior group probabilities}
  \item{assign}{a factor giving posterior group assignment}
  \item{tab}{matrix of retained principal components of PCA}
  \item{loadings}{principal axes of DAPC, giving coefficients of the linear
    combination of retained PCA axes.}
  \item{ind.coord}{principal components of DAPC, giving the coordinates
    of individuals onto principal axes of DAPC; also called the
    discriminant functions.}
  \item{grp.coord}{coordinates of the groups onto the principal axes of DAPC.}
  \item{posterior}{a data.frame giving posterior membership probabilities for
    all individuals and all clusters.}
  \item{var.contr}{(optional) a data.frame giving the contributions of original
    variables (alleles in the case of genetic data) to the principal components
    of DAPC.}
  \item{match.prp}{a list, where each item is the proportion of individuals 
    correctly matched to their original population in cross-validation.}

  
  === other outputs ===\cr
  Other functions have different outputs:\cr
  - \code{summary.dapc} returns a list with 6 components: \code{n.dim} (number
  of retained DAPC axes), \code{n.pop} (number of groups/populations),
  \code{assign.prop} (proportion of overall correct assignment),
  \code{assign.per.pop} (proportion of correct assignment per group),
  \code{prior.grp.size} (prior group sizes), and \code{post.grp.size} (posterior
  group sizes),  \code{xval.dapc}, \code{xval.genind} and \code{xval}
  (all return a list of four lists, each one with as many items as
  cross-validation runs.  The first item is a list of \code{assign}
  components, the secon is a list of \code{posterior} components, the
  thirs is a list of \code{ind.score} components and the fourth is a
  list of \code{match.prp} items, i.e. the prortion of the validation
  set correctly matched to its original population)
}
\references{
  Jombart T, Devillard S and Balloux F  (2010) Discriminant analysis of
  principal components: a new method for the analysis of genetically
  structured populations. BMC Genetics11:94. doi:10.1186/1471-2156-11-94
}
\seealso{
  \itemize{
  \item \code{\link{xvalDapc}}: selection of the optimal numbers of PCA axes
  retained in DAPC using cross-validation.

  \item \code{\link{scatter.dapc}}, \code{\link{assignplot}},
  \code{\link{compoplot}}: graphics for DAPC.

  \item \code{\link{find.clusters}}: to identify clusters without prior.

  \item \code{\link{dapcIllus}}: a set of simulated data illustrating
  the DAPC

  \item \code{\link{eHGDP}}, \code{\link{H3N2}}: empirical datasets
  illustrating DAPC
  }
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
## data(dapcIllus), data(eHGDP), and data(H3N2) illustrate the dapc
## see ?dapcIllus, ?eHGDP, ?H3N2
##
\dontrun{
example(dapcIllus)
example(eHGDP)
example(H3N2)
}

## H3N2 EXAMPLE ##
data(H3N2)
pop(H3N2) <- factor(H3N2$other$epid)
dapc1 <- dapc(H3N2, var.contrib=FALSE, scale=FALSE, n.pca=150, n.da=5)

## remove internal segments and ellipses, different pch, add MStree
scatter(dapc1, cell=0, pch=18:23, cstar=0, mstree=TRUE, lwd=2, lty=2)

## only ellipse, custom labels
scatter(dapc1, cell=2, pch="", cstar=0, posi.da="top",
lab=paste("year\n",2001:2006), axesel=FALSE, col=terrain.colors(10))


## SHOW COMPOPLOT ON MICROBOV DATA ##
data(microbov)
dapc1 <- dapc(microbov, n.pca=20, n.da=15)
compoplot(dapc1, lab="")




\dontrun{
## EXAMPLE USING GENLIGHT OBJECTS ##
## simulate data
x <- glSim(50,4e3-50, 50, ploidy=2)
x
plot(x)

## perform DAPC
dapc1 <- dapc(x, n.pca=10, n.da=1)
dapc1

## plot results
scatter(dapc1, scree.da=FALSE)

## SNP contributions
loadingplot(dapc1$var.contr)
loadingplot(tail(dapc1$var.contr, 100), main="Loading plot - last 100 SNPs")



## USE "PREDICT" TO PREDICT GROUPS OF NEW INDIVIDUALS ##
## load data
data(sim2pop)

## we make a dataset of:
## 30 individuals from pop A
## 30 individuals from pop B
## 30 hybrids

## separate populations and make F1
temp <- seppop(sim2pop)
temp <- lapply(temp, function(e) hybridize(e,e,n=30)) # force equal popsizes

## make hybrids
hyb <- hybridize(temp[[1]], temp[[2]], n=30)

## repool data - needed to ensure allele matching
newdat <- repool(temp[[1]], temp[[2]], hyb)
pop(newdat) <- rep(c("pop A", "popB", "hyb AB"), c(30,30,30))

## perform the DAPC on the first 2 pop (60 first indiv)
dapc1 <- dapc(newdat[1:60],n.pca=5,n.da=1)

## plot results
scatter(dapc1, scree.da=FALSE)

## make prediction for the 30 hybrids
hyb.pred <- predict(dapc1, newdat[61:90])
hyb.pred

## plot the inferred coordinates (circles are hybrids)
points(hyb.pred$ind.scores, rep(.1, 30))

## look at assignment using assignplot
assignplot(dapc1, new.pred=hyb.pred)
title("30 indiv popA, 30 indiv pop B, 30 hybrids")

## image using compoplot
compoplot(dapc1, new.pred=hyb.pred, ncol=2)
title("30 indiv popA, 30 indiv pop B, 30 hybrids")

## CROSS-VALIDATION ##
data(sim2pop)
xval <- xvalDapc(sim2pop@tab, pop(sim2pop), n.pca.max=100, n.rep=3)
xval
boxplot(xval$success~xval$n.pca, xlab="Number of PCA components",
ylab="Classification succes", main="DAPC - cross-validation")

}


}
\keyword{multivariate}
