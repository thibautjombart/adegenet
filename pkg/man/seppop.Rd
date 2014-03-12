\encoding{UTF-8}
\docType{methods}
\name{seppop}
\alias{seppop}
\alias{seppop-methods}
\alias{seppop,ANY-method}
\alias{seppop,genind-method}
\alias{seppop,genlight-method}
\title{ Separate genotypes per population}
\description{
  The function \code{seppop} splits a \linkS4class{genind} or a
  \linkS4class{genlight} object by population, returning a list of
  objects whose components each correspond to a population.\cr

  For \linkS4class{genind} objects, the output can either be a list of
  \linkS4class{genind} (default), or a list of matrices corresponding to
  the \code{@tab} slot.
}
\usage{
\S4method{seppop}{genind}(x,pop=NULL,truenames=TRUE,res.type=c("genind","matrix"),
  drop=FALSE, treatOther=TRUE, quiet=TRUE)

\S4method{seppop}{genlight}(x,pop=NULL, treatOther=TRUE, quiet=TRUE, \dots)
}
\arguments{
  \item{x}{a \linkS4class{genind} object}
  \item{pop}{a factor giving the population of each genotype in 'x'. If
    not provided, seeked in x\$pop.}
  \item{truenames}{a logical indicating whether true names should be
    used (TRUE, default) instead of generic labels (FALSE); used if
    res.type is "matrix".}
  \item{res.type}{a character indicating the type of returned results,
    a list of \linkS4class{genind} object (default) or a matrix of data
    corresponding to the 'tab' slots.}
  \item{drop}{a logical stating whether alleles that are no longer
    present in a subset of data should be discarded (TRUE) or kept anyway
    (FALSE, default).}
  \item{treatOther}{a logical stating whether elements of the
    \code{@other} slot should be treated as well (TRUE), or not
    (FALSE). See details in accessor documentations
    (\code{\link{pop}}).}
  \item{quiet}{a logical indicating whether warnings should be issued
    when trying to subset components of the \code{@other} slot (TRUE),
    or not (FALSE, default). }
  \item{\dots}{further arguments passed to the genlight constructor.}
}
\value{
  According to 'res.type': a list of \linkS4class{genind} object
  (default) or a matrix of data corresponding to the 'tab' slots.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\seealso{\code{\link{seploc}}, \code{\link{repool}}}
\examples{
\dontrun{
data(microbov)

obj <- seppop(microbov)
names(obj)

obj$Salers


#### example for genlight objects ####
x <- new("genlight", list(a=rep(1,1e3),b=rep(0,1e3),c=rep(1, 1e3)))
x

pop(x) # no population info
pop(x) <- c("pop1","pop2", "pop1") # set population memberships
pop(x)
seppop(x)
as.matrix(seppop(x)$pop1)[,1:20]
as.matrix(seppop(x)$pop2)[,1:20,drop=FALSE]
}
}
\keyword{manip}