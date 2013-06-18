\encoding{UTF-8}
\name{glPlot}
\alias{glPlot}
\alias{plot.genlight}
\alias{plot,genlight-method}
\title{Plotting genlight objects}
\description{
  \linkS4class{genlight} object can be plotted using the function
  \code{glPlot}, which is also used as the dedicated \code{plot}
  method. These functions relie on \code{\link{image}} to represent SNPs
  data. More specifically, colors are used to represent the number of
  second allele for each locus and individual.
}
\usage{
glPlot(x, col=NULL, legend=TRUE, posi="bottomleft", bg=rgb(1,1,1,.5),\dots)

\S4method{plot}{genlight}(x, y=NULL, col=NULL, legend=TRUE, posi="bottomleft", bg=rgb(1,1,1,.5),\dots)

}
\arguments{
  \item{x}{a \linkS4class{genlight} object.}
  \item{col}{an optional color vector; the first value corresponds to 0
    alleles, the last value corresponds to the ploidy level of the
    data. Therefore, the vector should have a length of (\code{ploidy(x)+1}).}
  \item{legend}{a logical indicating whether a legend should be added to
    the plot.}
  \item{posi}{a character string indicating where the legend should be
    positioned. Can be any concatenation of "bottom"/"top" and "left"/"right".}
  \item{bg}{a color used as a background for the legend; by default,
    transparent white is used; this may not be supported on some
    devices, and therefore background should be specified
    (e.g. \code{bg="white"}).}
  \item{\dots}{further arguments to be passed to \code{\link{image}}.}
  \item{y}{ununsed argument, present for compatibility with the
    \code{plot} generic.}
}
\seealso{
  - \code{\link{genlight}}: class of object for storing massive binary
  SNP data.
  
  - \code{\link{glSim}}: a simple simulator for \linkS4class{genlight}
  objects.

  - \code{\link{glPca}}: PCA for \linkS4class{genlight} objects.

}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
## simulate data
x <- glSim(100, 1e3, n.snp.struc=100, ploid=2)

## default plot
glPlot(x)
plot(x) # identical plot

## disable legend
plot(x, leg=FALSE)

## use other colors
plot(x, col=heat.colors(3), bg="white")
}
}
\keyword{multivariate}