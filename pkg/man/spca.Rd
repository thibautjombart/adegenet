\encoding{UTF-8}
\name{spca}
\alias{spca}
\alias{print.spca}
\alias{summary.spca}
\alias{plot.spca}
\alias{screeplot.spca}
\alias{colorplot.spca}
\title{Spatial principal component analysis}
\description{These functions are designed to perform a spatial principal
  component analysis and to display the results. They call upon
  \code{multispati} from the \code{ade4} package.
  
  \code{spca} performs the spatial component analysis. Other
  functions are:\cr
  
  - \code{print.spca}: prints the spca content\cr
  
  - \code{summary.spca}: gives variance and autocorrelation\cr
  statistics
  
  - \code{plot.spca}: usefull graphics (connection network, 3 different
  representations of map of scores, eigenvalues barplot and
  decomposition)\cr
  
  - \code{screeplot.spca}: decomposes spca eigenvalues into variance and
  autocorrelation\cr

  - \code{colorplot.spca}: represents principal components of sPCA in
  space using the RGB system.\cr

  A tutorial describes how to perform a sPCA: see
  \url{http://adegenet.r-forge.r-project.org/files/tutorial-spca.pdf} or
  type \code{adegenetTutorial(which="spca")}.
}
\usage{
spca(obj, xy=NULL, cn=NULL, matWeight=NULL,
     scale=FALSE, scannf=TRUE, nfposi=1, nfnega=1,
     type=NULL, ask=TRUE, plot.nb=TRUE, edit.nb=FALSE,
     truenames=TRUE, d1=NULL, d2=NULL, k=NULL, a=NULL, dmin=NULL)

\method{print}{spca}(x, \dots)

\method{summary}{spca}(object, \dots, printres=TRUE)

\method{plot}{spca}(x, axis = 1, useLag=FALSE, \dots)

\method{screeplot}{spca}(x, \dots, main=NULL)

\method{colorplot}{spca}(x, axes=1:ncol(x$li), useLag=FALSE, \dots)
}
\arguments{
  \item{obj}{a \code{genind} or \code{genpop} object.}
  \item{xy}{a matrix or data.frame with two columns for x and y
    coordinates. Seeked from obj\$other\$xy if it exists when xy is not
    provided. Can be NULL if a \code{nb} object is provided in
    \code{cn}.\cr
    Longitude/latitude coordinates should be converted first by a given
    projection (see 'See Also' section).}
  \item{cn}{a connection network of the class 'nb' (package spdep). Can
    be NULL if xy is provided. Can be easily obtained using the function
    chooseCN (see details).}
  \item{matWeight}{a square matrix of spatial weights, indicating the
    spatial proximities between entities. If provided, this argument
    prevails over \code{cn} (see details).}
  \item{scale}{a logical indicating whether alleles should be scaled to
    unit variance (TRUE) or not (FALSE, default).} 
  \item{scannf}{a logical stating whether eigenvalues should be chosen
    interactively (TRUE, default) or not (FALSE).}
  \item{nfposi}{an integer giving the number of positive eigenvalues retained
    ('global structures').}
  \item{nfnega}{an integer giving the number of negative eigenvalues retained
    ('local structures').}
  \item{type}{an integer giving the type of graph (see details in
    \code{chooseCN} help page). If provided, \code{ask} is set to FALSE.}
  \item{ask}{a logical stating whether graph should be chosen
    interactively (TRUE,default) or not (FALSE).}
  \item{plot.nb}{a logical stating whether the resulting graph should be
    plotted (TRUE, default) or not  (FALSE).}
  \item{edit.nb}{a logical stating whether the resulting graph should be
    edited manually for corrections (TRUE) or not  (FALSE, default).}  
  \item{truenames}{a logical stating whether true names should be used
    for 'obj' (TRUE, default) instead of generic labels (FALSE)}
  \item{d1}{the minimum distance between any two neighbours. Used if
    \code{type=5.}}
  \item{d2}{the maximum distance between any two neighbours. Used if
    \code{type=5}.}
  \item{k}{the number of neighbours per point. Used if
    \code{type=6}.}
  \item{a}{the exponent of the inverse distance matrix. Used if
    \code{type=7}.}
  \item{dmin}{the minimum distance between any two distinct points. Used
    to avoid infinite spatial proximities (defined as the inversed
    spatial distances). Used if \code{type=7}.}
  \item{x}{a \code{spca} object.}
  \item{object}{a \code{spca} object.}
  \item{printres}{a logical stating whether results should be printed on
    the screen (TRUE, default) or not (FALSE).}
  \item{axis}{an integer between 1 and (nfposi+nfnega) indicating which
    axis should be plotted.}
  \item{main}{a title for the screeplot; if NULL, a default one is
    used.}
  \item{\dots}{further arguments passed to other methods.}
  \item{axes}{the index of the columns of X to be represented. Up to
    three axes can be chosen.}
  \item{useLag}{a logical stating whether the lagged components
    (\code{x\$ls}) should be used instead of the components (\code{x\$li}).}
}
\details{
  The spatial principal component analysis (sPCA) is designed to
  investigate spatial patterns in the genetic variability. Given
  multilocus genotypes (individual level) or allelic frequency
  (population level) and spatial coordinates, it finds individuals (or
  population) scores maximizing the product of variance and spatial
  autocorrelation (Moran's I). Large positive and negative eigenvalues
  correspond to global and local structures.\cr

  Spatial weights can be obtained in several ways, depending how the
  arguments \code{xy}, \code{cn}, and \code{matWeight} are set.\cr
  When several acceptable ways are used at the same time, priority is as
  follows:\cr
  \code{matWeight} >  \code{cn} > \code{xy} \cr
}
\value{The class \code{spca} are given to lists with the following
  components:\cr
  \item{eig}{a numeric vector of eigenvalues.}
  \item{nfposi}{an integer giving the number of global structures
    retained.}
  \item{nfnega}{an integer giving the number of local structures retained.}
  \item{c1}{a data.frame of alleles loadings for each axis.}
  \item{li}{a data.frame of row (individuals or populations) coordinates
    onto the sPCA axes.}
  \item{ls}{a data.frame of lag vectors of the row coordinates; useful
    to clarify maps of global scores .}
  \item{as}{a data.frame giving the coordinates of the PCA axes onto the
    sPCA axes.}
  \item{call}{the matched call.}
  \item{xy}{a matrix of spatial coordinates.}
  \item{lw}{a list of spatial weights of class \code{listw}.}
  
  Other functions have different outputs:\cr
  - \code{summary.spca} returns a list with 3 components: \code{Istat}
  giving the null, minimum and maximum Moran's I values; \code{pca}
  gives variance and I statistics for the principal component analysis;
  \code{spca} gives variance and I statistics for the sPCA.\cr

  - \code{plot.spca} returns the matched call.\cr

  - \code{screeplot.spca} returns the matched call.  
}
\references{
Jombart, T., Devillard, S., Dufour, A.-B. and Pontier, D.
Revealing cryptic spatial patterns in genetic variability by a new
multivariate method. \emph{Heredity}, \bold{101}, 92--103.

Wartenberg, D. E. (1985) Multivariate spatial correlation: a method for
exploratory geographical analysis. \emph{Geographical Analysis},
\bold{17}, 263--283.

Moran, P.A.P. (1948) The interpretation of statistical
maps. \emph{Journal of the Royal Statistical Society, B}
\bold{10}, 243--251.

Moran, P.A.P. (1950) Notes on continuous stochastic
phenomena. \emph{Biometrika}, \bold{37}, 17--23.

de Jong, P. and Sprenger, C. and van Veen, F. (1984) On extreme values
of Moran's I and Geary's c. \emph{Geographical Analysis}, \bold{16}, 17--24.

}
\seealso{\code{\link{spcaIllus}} and \code{\link{rupica}} for datasets illustrating the sPCA \cr
  \code{\link{global.rtest}} and \code{\link{local.rtest}} \cr
  \code{\link{chooseCN}}, \code{\link[ade4]{multispati}},
  \code{\link[ade4]{multispati.randtest}}\cr
  \code{convUL}, from the package 'PBSmapping' to convert longitude/latitude to
  UTM coordinates.
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
## data(spcaIllus) illustrates the sPCA
## see ?spcaIllus
##
\dontrun{
example(spcaIllus)
example(rupica)
}
}
\keyword{spatial}
\keyword{multivariate}