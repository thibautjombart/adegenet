\name{colorplot}
\alias{colorplot}
\alias{colorplot.default}
\title{Represents a cloud of points with colors}
\description{
  The \code{colorplot} function represents a cloud of points with colors
  corresponding to a combination of 1,2 or 3 quantitative variables,
  assigned to RGB (Red, Green, Blue) channels. For instance, this can be useful to
  represent up to 3 principal components in space. Note that the
  property of such representation to convey multidimensional information
  has not been investigated.\cr

  \code{colorplot} is a S3 generic function. Methods are defined for
  particular objects, like \code{\link{spca}} objects.

}
\usage{
colorplot(\dots)

\method{colorplot}{default}(xy, X, axes=NULL, add.plot=FALSE, defaultLevel=0, transp=FALSE, alpha=.5, \dots)
}
\arguments{
  \item{xy}{a numeric matrix with two columns (e.g. a matrix of spatial coordinates.}
  \item{X}{a matrix-like containing numeric values that are translated
    into the RGB system. Variables are considered to be in columns.}
  \item{axes}{the index of the columns of X to be represented. Up to
    three axes can be chosen. If null, up to the first three columns of
    X are used.}
  \item{add.plot}{a logical stating whether the colorplot should be
    added to the existing plot (defaults to FALSE).}
  \item{defaultLevel}{a numeric value between 0 and 1, giving the
    default level in a color for which values are not specified. Used
    whenever less than three axes are specified.}
  \item{transp}{a logical stating whether the produced colors should be
    transparent (TRUE) or not (FALSE, default).}
  \item{alpha}{the alpha level for transparency, between 0 (fully
    transparent) and 1 (not transparent); see \code{?rgb} for more details.}
  \item{\dots}{further arguments to be passed to other methods. In
    \code{colorplot.default}, these arguments are passed to plot/points
    functions. See \code{?plot.default} and \code{?points}.}
}
\value{
  Invisibly returns a vector of colours used in the plot.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\examples{
# a toy example
xy <- expand.grid(1:10,1:10)
df <- data.frame(x=1:100, y=100:1, z=runif(100,0,100))
colorplot(xy,df,cex=10,main="colorplot: toy example")

\dontrun{
# a genetic example using a sPCA
if(require(spdep)){
data(spcaIllus)
dat3 <- spcaIllus$dat3
spca3 <- spca(dat3,xy=dat3$other$xy,ask=FALSE,type=1,plot=FALSE,scannf=FALSE,nfposi=1,nfnega=1)
colorplot(spca3, cex=4, main="colorplot: a sPCA example")
text(spca3$xy[,1], spca3$xy[,2], dat3$pop)
mtext("P1-P2 in cline\tP3 random \tP4 local repulsion")
}
}
}
\keyword{multivariate}
\keyword{hplot}
