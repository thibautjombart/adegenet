\name{loadingplot}
\alias{loadingplot}
\alias{loadingplot.default}
\title{Represents a cloud of points with colors}
\description{
  The \code{loadingplot} function represents positive values of a vector
  and identifies the values above a given threshold. It can also
  indicate groups of observations provided as a factor. \cr

  Such graphics can be used, for instance, to assess the weight of each
  variable (loadings) in a given analysis.
}
\usage{
loadingplot(x, \dots)

\method{loadingplot}{default}(x, at=NULL, threshold=quantile(x,0.75),
            axis=1, fac=NULL, byfac=FALSE,
            lab=NULL, cex.lab=0.7, cex.fac=1, lab.jitter=0,
            main="Loading plot", xlab="Variables", ylab="Loadings",
            srt = 0, adj = NULL, \dots)
}
\arguments{
  \item{x}{either a vector with numeric values to be plotted, or a
    matrix-like object containing numeric values. In such case, the
    \code{x[,axis]} is used as vector of values to be plotted.}
  \item{at}{an optional numeric vector giving the abscissa at which
    loadings are plotted. Useful when variates are SNPs with a known
    position in an alignement.}
  \item{threshold}{a threshold value above which values of x are
    identified. By default, this is the third quartile of x.}
  \item{axis}{an integer indicating the column of x to be plotted; used
    only if x is a matrix-like object.}
  \item{fac}{a factor defining groups of observations.}
  \item{byfac}{a logical stating whether loadings should be averaged by
    groups of observations, as defined by \code{fac}.}
  \item{lab}{a character vector giving the labels used to annotate
    values above the threshold; if NULL, names are taken from the object.}
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
  \item{\dots}{further arguments to be passed to the plot function.}
}
  \value{
    Invisibly returns a list with the following components:\cr
    - threshold: the threshold used\cr
    - var.names: the names of observations above the threshold\cr
    - var.idx: the indices of observations above the threshold\cr
    - var.values: the values above the threshold\cr
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\examples{
x <- runif(20)
names(x) <- letters[1:20]
grp <- factor(paste("group", rep(1:4,each=5)))

## basic plot
loadingplot(x)

## adding groups
loadingplot(x,fac=grp,main="My title",cex.lab=1)
}
\keyword{multivariate}
\keyword{hplot}
