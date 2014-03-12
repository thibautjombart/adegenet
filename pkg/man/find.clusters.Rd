\encoding{UTF-8}
\name{find.clusters}
\alias{find.clusters}
\alias{find.clusters.data.frame}
\alias{find.clusters.matrix}
\alias{find.clusters.genind}
\alias{find.clusters.genlight}
\alias{.find.sub.clusters}
\title{find.cluster: cluster identification using successive K-means}
\description{
   These functions implement the clustering procedure used in
  Discriminant Analysis of Principal Components (DAPC, Jombart et
  al. 2010). This procedure consists in running successive K-means with
  an increasing number of clusters (\code{k}), after transforming data
  using a principal component analysis (PCA). For each model, a
  statistical measure of goodness of fit (by default, BIC) is computed,
  which allows to choose the optimal \code{k}. See \code{details} for a
  description of how to select the optimal \code{k} and
  \code{vignette("adegenet-dapc")} for a tutorial.

  Optionally, hierarchical clustering can be sought by providing a prior
  clustering of individuals (argument \code{clust}). In such case, clusters will
  be sought within each prior group.

  The K-means procedure used in \code{find.clusters} is \code{\link[stats]{kmeans}} function
  from the \code{stats} package. The PCA function is \code{\link[ade4]{dudi.pca}} from the
  \code{ade4} package, except for \linkS4class{genlight} objects which
  use the \code{\link{glPca}} procedure from adegenet.

  \code{find.clusters} is a generic function with methods for the
 following types of objects:\cr
 - \code{data.frame} (only numeric data)\cr
 - \code{matrix} (only numeric data)\cr
 - \code{\linkS4class{genind}} objects (genetic markers)\cr
 - \code{\linkS4class{genlight}} objects (genome-wide SNPs)
 
}
\usage{
\method{find.clusters}{data.frame}(x, clust=NULL, n.pca=NULL,
              n.clust=NULL, stat=c("BIC","AIC", "WSS"),
              choose.n.clust=TRUE,criterion=c("diffNgroup",
              "min","goesup", "smoothNgoesup", "goodfit"),
              max.n.clust=round(nrow(x)/10), n.iter=1e5, n.start=10,
              center=TRUE, scale=TRUE, pca.select=c("nbEig","percVar"),
              perc.pca=NULL, \ldots, dudi=NULL)

\method{find.clusters}{matrix}(x, \ldots)

\method{find.clusters}{genind}(x, clust=NULL, n.pca=NULL, n.clust=NULL,
              stat=c("BIC","AIC", "WSS"), choose.n.clust=TRUE,
              criterion=c("diffNgroup", "min","goesup", "smoothNgoesup",
              "goodfit"), max.n.clust=round(nrow(x@tab)/10), n.iter=1e5,
              n.start=10, scale=FALSE, truenames=TRUE, \ldots)

\method{find.clusters}{genlight}(x, clust=NULL, n.pca=NULL,
              n.clust=NULL, stat=c("BIC", "AIC",
              "WSS"),choose.n.clust=TRUE, criterion=c("diffNgroup",
              "min","goesup","smoothNgoesup", "goodfit"),
              max.n.clust=round(nInd(x)/10), n.iter=1e5,n.start=10,
              scale=FALSE, pca.select=c("nbEig","percVar"),
              perc.pca=NULL,glPca=NULL, \ldots)
}
\arguments{
  \item{x}{\code{a data.frame}, \code{matrix}, or \code{\linkS4class{genind}}
    object. For the \code{data.frame} and \code{matrix} arguments, only
    quantitative variables should be provided.}
  \item{clust}{an optional \code{factor} indicating a prior group membership of
    individuals. If provided, sub-clusters will be sought within each prior
    group.}
  \item{n.pca}{an \code{integer} indicating the number of axes retained in the
    Principal Component Analysis (PCA) step. If \code{NULL}, interactive selection
    is triggered.}
  \item{n.clust}{ an optinal \code{integer} indicating the number of clusters to
    be sought. If provided, the function will only run K-means once, for this
    number of clusters. If left as \code{NULL}, several K-means are run for a
    range of k (number of clusters) values.}
  \item{stat}{ a \code{character} string matching 'BIC', 'AIC', or 'WSS', which
    indicates the statistic to be computed for each model (i.e., for each value of
    \code{k}). BIC: Bayesian Information Criterion. AIC: Aikaike's Information
    Criterion. WSS: within-groups sum of squares, that is, residual variance.}
  \item{choose.n.clust}{ a \code{logical} indicating whether the number of
    clusters should be chosen by the user (TRUE, default), or automatically,
    based on a given criterion (argument \code{criterion}). It is HIGHLY
    RECOMMENDED to choose the number of clusters INTERACTIVELY, since
    i) the decrease of the summary statistics (BIC by default) is
    informative, and ii) no criteria for automatic selection is
    appropriate to all cases (see details).}
  \item{criterion}{ a \code{character} string matching "diffNgroup",
    "min","goesup", "smoothNgoesup", or "conserv", indicating the criterion for automatic
    selection of the optimal number of clusters. See \code{details} for
    an explanation of these procedures.}
  \item{max.n.clust}{ an \code{integer} indicating the maximum number of
    clusters to be tried. Values of 'k' will be picked up between 1 and \code{max.n.clust}}
  \item{n.iter}{ an \code{integer} indicating the number of iterations to be used
    in each run of K-means algorithm. Corresponds to \code{iter.max} of
    \code{kmeans} function.}
  \item{n.start}{ an \code{integer} indicating the number of randomly
    chosen starting centroids to be used in each run of the K-means
    algorithm. Using more starting points ensures convergence of the
    algorithm. Corresponds to \code{nstart} of \code{kmeans} function.}
  \item{center}{a \code{logical} indicating whether variables should be centred to
    mean 0 (TRUE, default) or not (FALSE). Always TRUE for \linkS4class{genind}
    objects.}
  \item{scale}{a \code{logical} indicating whether variables should be
    scaled (TRUE) or not (FALSE, default). Scaling consists in dividing
    variables by their (estimated) standard deviation to account for
    trivial differences in variances. In allele frequencies, it comes with
    the risk of giving uninformative alleles more importance while
    downweighting informative alleles. Further scaling options are
    available for \linkS4class{genind} objects (see argument
    \code{scale.method}).}
  \item{pca.select}{a \code{character} indicating the mode of selection of PCA
    axes, matching either "nbEig" or "percVar". For "nbEig", the user
    has to specify the number of axes retained (interactively, or via
    \code{n.pca}). For "percVar", the user has to specify the minimum amount of
    the total variance to be preserved by the retained axes, expressed as a
    percentage (interactively, or via \code{perc.pca}).  }
  \item{perc.pca}{a \code{numeric} value between 0 and 100 indicating the
    minimal percentage of the total variance of the data to be expressed by the
    retained axes of PCA.} 
  \item{truenames}{a \code{logical} indicating whether true (i.e., user-specified)
    labels should be used in object outputs (TRUE, default) or not
    (FALSE), in which case generic labels are used.}
  \item{\ldots}{further arguments to be passed to other functions. For
    \code{find.clusters.matrix}, arguments are to match those of the
    \code{data.frame} method.}
  \item{dudi}{optionally, a multivariate analysis with the class
    \code{dudi} (from the ade4 package). If provided, prior PCA will be
    ignored, and this object will be used as a prior step for variable
    orthogonalisation.}
  \item{glPca}{an optional \code{\link{glPca}} object; if provided,
    dimension reduction is not performed (saving computational time) but
    taken directly from this object.}
}
\details{
  === ON THE SELECTION OF K ===\cr
  (where K is the 'optimal' number of clusters)

  So far, the analysis of data simulated under various population
  genetics models (see reference) suggested an ad hoc rule for the
  selection of the optimal number of clusters. First important result is
  that BIC seems for efficient than AIC and WSS to select the
  appropriate number of clusters (see example). The rule of thumb
  consists in increasing K until it no longer leads to an appreciable
  improvement of fit (i.e., to a decrease of BIC).  In the most simple
  models (island models), BIC decreases until it reaches the optimal K,
  and then increases. In these cases, our rule amounts to choosing the
  lowest K. In other models such as stepping stones, the decrease of BIC
  often continues after the optimal K, but is much less steep.


  An alternative approach is the automatic selection based on a fixed
  criterion. Note that, in any case, it is highly recommended to look at
  the graph of the BIC for different numbers of clusters as displayed
  during the interactive cluster selection.
  To use automated selection, set \code{choose.n.clust} to FALSE and specify
  the \code{criterion} you want to use, from the following values:

  - "diffNgroup": differences between successive values of the summary
  statistics (by default, BIC) are splitted into two groups using a
  Ward's clustering method (see \code{?hclust}), to differentiate sharp
  decrease from mild decreases or increases. The retained K is the one
  before the first group switch. Appears to work well for
  island/hierarchical models, and decently for isolation by distance
  models, albeit with some unstability. Can be impacted by an initial,
  very sharp decrease of the test statistics. IF UNSURE ABOUT THE
  CRITERION TO USE, USE THIS ONE.

  - "min": the model with the minimum summary statistics (as specified
  by \code{stat} argument, BIC by default) is retained. Is likely to
  work for simple island model, using BIC. It is likely to fail in
  models relating to stepping stones, where the BIC always decreases
  (albeit by a small amount) as K increases. In general, this approach
  tends to over-estimate the number of clusters.

  - "goesup": the selected model is the K after which increasing the
  number of clusters leads to increasing the summary statistics. Suffers
  from inaccuracy, since i) a steep decrease might follow a small 'bump'
  of increase of the statistics, and ii) increase might never happen, or
  happen after negligible decreases. Is likely to work only for
  clear-cut island models.

  - "smoothNgoesup": a variant of "goesup", in which the summary
  statistics is first smoothed using a lowess approach. Is meant to be
  more accurate than "goesup" as it is less prone to stopping to small
  'bumps' in the decrease of the statistics.

  - "goodfit": another criterion seeking a good fit with a minimum
  number of clusters. This approach does not rely on differences between
  successive statistics, but on absolute fit. It selects the model with
  the smallest K so that the overall fit is above a given threshold.
}
\value{
  The class \code{find.clusters} is a list with the following
  components:\cr
  \item{Kstat}{a \code{numeric} vector giving the values of the summary
  statistics for the different values of K. Is NULL if \code{n.clust} was specified.}
  \item{stat}{a \code{numeric} value giving the value of the summary statistics
  for the retained model}
  \item{grp}{a \code{factor} giving group membership for each individual.}
  \item{size}{an \code{integer} vector giving the size of the different clusters.}
}
\references{
  Jombart T, Devillard S and Balloux F (2010) Discriminant analysis of
  principal components: a new method for the analysis of genetically
  structured populations. BMC Genetics 11:94. doi:10.1186/1471-2156-11-94
}
\seealso{
  - \code{\link{dapc}}: implements the DAPC.

  - \code{\link{scatter.dapc}}: graphics for DAPC.

  - \code{\link{dapcIllus}}: dataset illustrating the DAPC and \code{find.clusters}.

  - \code{\link{eHGDP}}: dataset illustrating the DAPC and \code{find.clusters}.

  - \code{\link[stats]{kmeans}}: implementation of K-means in the stat package.

  - \code{\link[ade4]{dudi.pca}}: implementation of PCA in the ade4 package.

}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
## THIS ONE TAKES A FEW MINUTES TO RUN ## 
data(eHGDP)

## here, n.clust is specified, so that only on K value is used
grp <- find.clusters(eHGDP, max.n=30, n.pca=200, scale=FALSE,
n.clust=4) # takes about 2 minutes
names(grp)
grp$Kstat
grp$stat


## to try different values of k (interactive)
grp <- find.clusters(eHGDP, max.n=50, n.pca=200, scale=FALSE)

## and then, to plot BIC values:
plot(grp$Kstat, type="b", col="blue")



## ANOTHER SIMPLE EXAMPLE ## 
data(sim2pop) # this actually contains 2 pop

## DETECTION WITH BIC (clear result)
foo.BIC <- find.clusters(sim2pop, n.pca=100, choose=FALSE)
plot(foo.BIC$Kstat, type="o", xlab="number of clusters (K)", ylab="BIC",
col="blue", main="Detection based on BIC")
points(2, foo.BIC$Kstat[2], pch="x", cex=3)
mtext(3, tex="'X' indicates the actual number of clusters")


## DETECTION WITH AIC (less clear-cut)
foo.AIC <- find.clusters(sim2pop, n.pca=100, choose=FALSE, stat="AIC")
plot(foo.AIC$Kstat, type="o", xlab="number of clusters (K)",
ylab="AIC", col="purple", main="Detection based on AIC")
points(2, foo.AIC$Kstat[2], pch="x", cex=3)
mtext(3, tex="'X' indicates the actual number of clusters")


## DETECTION WITH WSS (less clear-cut)
foo.WSS <- find.clusters(sim2pop, n.pca=100, choose=FALSE, stat="WSS")
plot(foo.WSS$Kstat, type="o", xlab="number of clusters (K)", ylab="WSS
(residual variance)", col="red", main="Detection based on WSS")
points(2, foo.WSS$Kstat[2], pch="x", cex=3)
mtext(3, tex="'X' indicates the actual number of clusters")


## TOY EXAMPLE FOR GENLIGHT OBJECTS ##
x <- glSim(100,500,500)
x
plot(x)
grp <- find.clusters(x, n.pca=100, choose=FALSE, stat="BIC")
plot(grp$Kstat, type="o", xlab="number of clusters (K)",
ylab="BIC",main="find.clusters on a genlight object\n(two groups)")
}
}
\keyword{multivariate}
