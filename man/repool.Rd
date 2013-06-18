\name{repool}
\alias{repool}
\title{Pool several genotypes into a single dataset}
\description{
  The function \code{repool} allows to merge genotypes from different
  \linkS4class{genind} objects into a single 'pool' (i.e. a new \linkS4class{genind}).
  The markers have to be the same for all objects to be merged, but
  there is no constraint on alleles.\cr

  This function can be useful, for instance, when hybrids are created
  using \code{\link{hybridize}}, to merge hybrids with their parent
  population for further analyses. Note that \code{repool} can also
  reverse the action of \code{\link{seppop}}.
}
\usage{
repool(\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{can be i) a list whose components are valid
    \linkS4class{genind} objects or, ii) several valid
    \linkS4class{genind} objects separated by commas.}
}
\value{
  A \linkS4class{genind} object.
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\seealso{\code{\link{seploc}}, \code{\link{seppop}}}
\examples{
\dontrun{
## use the cattle breeds dataset
data(microbov)
temp <- seppop(microbov)
names(temp)

## hybridize salers and zebu -- nasty cattle
zebler <- hybridize(temp$Salers, temp$Zebu, n=40)
zebler

## now merge zebler with other cattle breeds
nastyCattle <- repool(microbov, zebler)
nastyCattle
}
}
\keyword{manip}
