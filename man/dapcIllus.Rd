\encoding{UTF-8}
\name{dapcIllus}
\alias{dapcIllus}
\docType{data}
\title{Simulated data illustrating the DAPC}
\description{
  Datasets illustrating the Discriminant Analysis of Principal
  Components (DAPC, Jombart et al. submitted).\cr

  These data were simulated using various models using Easypop (2.0.1).
  The \code{dapcIllus} is a list containing the
  following \linkS4class{genind} objects:\cr
  - "a": island model with 6 populations \cr
  - "b": hierarchical island model with 6 populations (3,2,1) \cr
  - "c": one-dimensional stepping stone with 2x6 populations, and a boundary between the
  two sets of 6 populations\cr
  - "d": one-dimensional stepping stone with 24 populations\cr

  See "source" for a reference providing simulation details.
}
\usage{data(dapcIllus)}
\format{
  \code{dapcIllus} is list of 4 components being all genind objects.
}
\source{
  Jombart, T., Devillard, S. and Balloux, F.  Discriminant analysis of
  principal components: a new method for the analysis of genetically
  structured populations. Submitted to \emph{BMC genetics}.

}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}
\seealso{
  - \code{\link{dapc}}: implements the DAPC.

  - \code{\link{eHGDP}}: dataset illustrating the DAPC and \code{find.clusters}.

  - \code{\link{H3N2}}: dataset illustrating the DAPC.

  - \code{\link{find.clusters}}: to identify clusters without prior.
}
\references{
   Jombart, T., Devillard, S. and Balloux, F.  Discriminant analysis of
   principal components: a new method for the analysis of genetically
   structured populations. Submitted to \emph{Genetics}.
}
\examples{
\dontrun{

data(dapcIllus)
attach(dapcIllus)
a # this is a genind object, like b, c, and d.


## FINS CLUSTERS EX NIHILO
clust.a <- find.clusters(a, n.pca=100, n.clust=6)
clust.b <- find.clusters(b, n.pca=100, n.clust=6)
clust.c <- find.clusters(c, n.pca=100, n.clust=12)
clust.d <- find.clusters(d, n.pca=100, n.clust=24)

## examin outputs
names(clust.a)
lapply(clust.a, head)


## PERFORM DAPCs
dapc.a <- dapc(a, pop=clust.a$grp, n.pca=100, n.da=5)
dapc.b <- dapc(b, pop=clust.b$grp, n.pca=100, n.da=5)
dapc.c <- dapc(c, pop=clust.c$grp, n.pca=100, n.da=11)
dapc.d <- dapc(d, pop=clust.d$grp, n.pca=100, n.da=23)


## LOOK AT ONE RESULT
dapc.a
summary(dapc.a)

## FORM A LIST OF RESULTS FOR THE 4 DATASETS
lres <- list(dapc.a, dapc.b, dapc.c, dapc.d)


## DRAW 4 SCATTERPLOTS
par(mfrow=c(2,2))
lapply(lres, scatter)


# detach data
detach(dapcIllus)
}
}
\keyword{datasets}
