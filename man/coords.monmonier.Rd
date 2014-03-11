\name{coords.monmonier}
\alias{coords.monmonier}
\title{Returns original points in results paths of an object of class 'monmonier'}
\description{
  The original implementation of \code{monmonier} in package
  \bold{adegenet} returns path coordinates, \code{coords.monmonier}
  additionally displays identities of the original points of the network,
  based on original coordinates.
}
\usage{
coords.monmonier(x)
}
\arguments{
  \item{x}{an object of class \code{\link{monmonier}}.}
}
\value{
  Returns a list with elements according to the \code{x$nrun} result of the \code{\link{monmonier}} object.
  Corresponding path points are in the same order as in the original object.

  run1 (run2, ...): for each run, a list containing a matrix giving the
  original points in the network (\code{first} and \code{second},
  indicating pairs of neighbours). Path coordinates are stored in
  columns \code{x.hw} and \code{y.hw}. \code{first} and \code{second}
  are integers referring to the row numbers in the \code{x$xy} matrix of
  the original \code{\link{monmonier}} object.
}
\author{Peter Solymos, \email{Solymos.Peter@aotk.szie.hu}, \url{http://www.univet.hu/users/psolymos/personal/}}
\seealso{
  \code{\link{monmonier}}
}
\examples{

\dontrun{
if(require(spdep)){

load(system.file("files/mondata1.rda",package="adegenet"))
cn1 <- chooseCN(mondata1$xy,type=2,ask=FALSE)
mon1 <- monmonier(mondata1$xy,dist(mondata1$x1),cn1,threshold=2,nrun=3)

mon1$run1
mon1$run2
mon1$run3
path.coords <- coords.monmonier(mon1)
path.coords
}
}
}
\keyword{methods}
\keyword{manip}
