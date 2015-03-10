\encoding{UTF-8}
\name{glPca}
\alias{glPca}
\alias{print.glPca}
\alias{scatter.glPca}
\alias{loadingplot.glPca}
\title{Principal Component Analysis for genlight objects}
\description{
  These functions implement Principal Component Analysis (PCA) for
  massive SNP datasets stored as \linkS4class{genlight} object. This
  implementation has the advantage of never representing to complete data
  matrix, therefore making huge economies in terms of rapid access
  memory (RAM). When the \code{parallel} package is available,
  \code{glPca} uses multiple-core ressources for more efficient
  computations. \code{glPca} returns lists with the class \code{glPca}
  (see 'value'). 
  
  Other functions are defined for objects of this class:
  
  - \code{print}: prints the content of a \code{glPca} object.
  
  - \code{scatter}: produces scatterplots of principal components,
  with a screeplot of eigenvalues as inset.
    
  - \code{loadingplot}: plots the loadings of the analysis for one given
  axis, using an adapted version of the generic function \code{loadingplot}.
}
\usage{
glPca(x, center = TRUE, scale = FALSE, nf = NULL, loadings = TRUE, 
    alleleAsUnit = FALSE, useC = TRUE, parallel = require("parallel"),
  n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)

\method{print}{glPca}(x, \dots)

\method{scatter}{glPca}(x, xax = 1, yax = 2, posi = "bottomleft", bg = "white", 
    ratio = 0.3, label = rownames(x$scores), clabel = 1, xlim = NULL, 
    ylim = NULL, grid = TRUE, addaxes = TRUE, origin = c(0, 0), 
    include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft", 
    cgrid = 1, pixmap = NULL, contour = NULL, area = NULL, \dots)

\method{loadingplot}{glPca}(x, at=NULL, threshold=NULL, axis=1,
    fac=NULL, byfac=FALSE, lab=rownames(x$loadings), cex.lab=0.7, cex.fac=1,
    lab.jitter=0, main="Loading plot", xlab="SNP positions",
    ylab="Contributions", srt = 90, adj = c(0, 0.5), \dots)

}
\arguments{
  \item{x}{for \code{glPca}, a \linkS4class{genlight} object; for
    \code{print}, \code{scatter}, and \code{loadingplot}, a \code{glPca}
  object.}
  \item{center}{a logical indicating whether the numbers of alleles
    should be centered; defaults to TRUE}
  \item{scale}{a logical indicating whether the numbers of alleles
    should be scaled; defaults to FALSE}
  \item{nf}{an integer indicating the number of principal components to
    be retained; if NULL, a screeplot of eigenvalues will be displayed
    and the user will be asked for a number of retained axes.}
  \item{loadings}{a logical indicating whether loadings of the alleles
    should be computed (TRUE, default), or not (FALSE). Vectors of
    loadings are not always useful, and can take a large amount of
    RAM when millions of SNPs are considered.}
  \item{alleleAsUnit}{a logical indicating whether alleles are
    considered as units (i.e., a diploid genotype equals two samples, a
    triploid, three, etc.) or whether individuals are considered as
    units of information.}
  \item{useC}{a logical indicating whether compiled C code should be
    used for faster computations; this option cannot be used alongside
    parallel option.}
  \item{parallel}{a logical indicating whether multiple cores -if
    available- should be used for the computations (TRUE, default), or
    not (FALSE); requires the package \code{parallel} to be installed
    (see details); this option cannot be used alongside useCoption.}
  \item{n.cores}{if \code{parallel} is TRUE, the number of cores to be
    used in the computations; if NULL, then the maximum number of cores
    available on the computer is used.}
  \item{returnDotProd}{a logical indicating whether the matrix of dot
    products between individuals should be returned (TRUE) or not
    (FALSE, default).}
  \item{matDotProd}{an optional matrix of dot products between
    individuals, NULL by default. This option is used internally to
    speed up computation time when re-running the same PCA several
    times. Leave this argument as NULL unless you really know what you
    are doing.}
  \item{\dots}{further arguments to be passed to other functions.}

  \item{xax,yax}{\code{integers} specifying which principal components
    should be shown in x and y axes. }
  \item{posi,bg,ratio}{arguments used to customize the inset in scatterplots
    of \code{glPca} results. See \code{\link[ade4]{add.scatter}} documentation in the
    ade4 package for more details.}
  \item{label,clabel,xlim,ylim,grid,addaxes,origin,include.origin,sub,csub,possub,cgrid,pixmap,contour,area}{
    arguments passed to \code{\link[ade4]{s.class}}; see \code{?s.label} for more information}
   \item{at}{an optional numeric vector giving the abscissa at which
    loadings are plotted. Useful when variates are SNPs with a known
    position in an alignement.}
  \item{threshold}{a threshold value above which values of x are
    identified. By default, this is the third quartile of x.}
  \item{axis}{an integer indicating the column of x to be plotted; used
    only if x is a matrix-like object.}
  \item{fac}{a factor defining groups of SNPs.}
  \item{byfac}{a logical stating whether loadings should be averaged by
    groups of SNPs, as defined by \code{fac}.}
  \item{lab}{a character vector giving the labels used to annotate
    values above the threshold.}
  \item{cex.lab}{a numeric value indicating the size of annotations.}
  \item{cex.fac}{a numeric value indicating the size of annotations for
    groups of observations.}
  \item{lab.jitter}{a numeric value indicating the factor of
    randomisation for the position of annotations. Set to 0 (by default)
  implies no randomisation.}
  \item{main}{the main title of the figure.}
  \item{xlab}{the title of the x axis.}
  \item{ylab}{the title of the y axis.}
  \item{srt}{rotation of the labels; see ?text.}
  \item{adj}{adjustment of the labels; see ?text.}
}
\details{
  === Using multiple cores ===
  
  Most recent machines have one or several processors with multiple
  cores. R processes usually use one single core. The package
  \code{parallel} allows for parallelizing some computations on
  multiple cores, which can decrease drastically computational time.

  Lastly, note that using compiled C code (\code{useC=TRUE})is an
  alternative for speeding up computations, but cannot be used together
  with the parallel option.
}
\value{
  === glPca objects ===
  
  The class \code{glPca} is a list with the following
  components:\cr
  \item{call}{the matched call.}
  \item{eig}{a numeric vector of eigenvalues.}
  \item{scores}{a matrix of principal components, containing the
    coordinates of each individual (in row) on each principal axis (in column).}
  \item{loadings}{(optional) a matrix of loadings, containing the loadings of each
    SNP (in row) for each principal axis (in column).}

  -
  
  === other outputs ===
  
  Other functions have different outputs:\cr
  - \code{scatter} return the matched call.\cr
  - \code{loadingplot} returns information about the most contributing
  SNPs (see \code{\link{loadingplot.default}})
}

\seealso{
  - \code{\linkS4class{genlight}}: class of object for storing massive binary
  SNP data.

  - \code{\link{glSim}}: a simple simulator for \linkS4class{genlight} objects.

  - \code{\link{glPlot}}: plotting \linkS4class{genlight} objects.

  - \code{\link{dapc}}: Discriminant Analysis of Principal Components.

}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
## simulate a toy dataset
x <- glSim(50,4e3, 50, ploidy=2)
x
plot(x)

## perform PCA
pca1 <- glPca(x, nf=2)

## plot eigenvalues
barplot(pca1$eig, main="eigenvalues", col=heat.colors(length(pca1$eig)))

## basic plot
scatter(pca1, ratio=.2)

## plot showing groups
s.class(pca1$scores, pop(x), col=colors()[c(131,134)])
add.scatter.eig(pca1$eig,2,1,2)
}
}
\keyword{multivariate}