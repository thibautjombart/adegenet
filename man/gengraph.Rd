\encoding{UTF-8}
\name{gengraph}
\alias{gengraph}
\alias{gengraph.default}
\alias{gengraph.matrix}
\alias{gengraph.dist}
\alias{gengraph.genind}
\alias{gengraph.genpop}
\alias{gengraph.DNAbin}
\title{Genetic transitive graphs}
\description{
  These functions are under development. Please email the author before
  using them for published work.\cr

  The function \code{gengraph} generates graphs based on genetic
  distances, so that pairs of entities (individuals or populations) are
  connected if and only if they are distant by less than a given
  threshold distance. Graph algorithms and classes from the
  \code{\link[igraph]{igraph}} package are used.\cr

  \code{gengraph} is a generic function with methods for the
 following types of objects:\cr
 - \code{matrix} (only numeric data)\cr
 - \code{dist} \cr
 - \code{\linkS4class{genind}} objects (genetic markers, individuals)\cr
 - \code{\linkS4class{genpop}} objects (genetic markers, populations)\cr
 -  \code{\link[ape]{DNAbin}} objects (DNA sequences)
}
\usage{
gengraph(x, \dots)

\method{gengraph}{matrix}(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE,
plot=TRUE, show.graph=TRUE, col.pal=funky, truenames=TRUE, \dots)

\method{gengraph}{dist}(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE,
plot=TRUE, show.graph=TRUE, col.pal=funky, truenames=TRUE, \dots)

\method{gengraph}{genind}(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE,
plot=TRUE, show.graph=TRUE, col.pal=funky, truenames=TRUE, \dots)

\method{gengraph}{genpop}(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE,
         plot=TRUE, show.graph=TRUE, col.pal=funky, method=1,
         truenames=TRUE, \dots)

\method{gengraph}{DNAbin}(x, cutoff=NULL, ngrp=NULL, computeAll=FALSE,
plot=TRUE, show.graph=TRUE, col.pal=funky, truenames=TRUE, \dots)

}
\arguments{
  \item{x}{a \code{matrix}, \code{dist}, \code{\linkS4class{genind}},
    \code{\linkS4class{genpop}}, or \code{DNAbin} object. For
    \code{matrix} and \code{dist}, the object represents pairwise
    (by default, Hamming) distances between considered individuals.}
  \item{cutoff}{a \code{numeric} value indicating the cutoff point,
    i.e. the distance at which two entities are no longer connected in
    the garph produced by the method.}
  \item{ngrp}{an \code{integer} indicating the number of groups to be
    looked for. A message is issued if this exact number could not be found.}
  \item{computeAll}{a \code{logical} stating whether to investigate
    solutions for every (integer) cutoff point; defaults to FALSE.}
  \item{plot}{a \code{logical} indicating whether plots should be drawn;
    defaults to TRUE; this operation can take time for large, highly-connected graphs.}
  \item{show.graph}{a \code{logical} indicating whether the found graph
    should be drawn, only used in the interactive mode; this operation
    can take time for large, highly-connected graphs; defaults to FALSE.}
  \item{col.pal}{a color palette used to define group colors.}
  \item{method}{an \code{integer} ranging from 1 to 6 indicating the
    type of method to be used to derive a matrix of pairwise distances
    between populations; values from 1 to 5 are passed to the function
    \code{dist.genpop}; 6 corresponds to pairwise Fst; other values are
    not supported.}
  \item{truenames}{a logical indicating whether original labels should
    be used for plotting (TRUE), as opposed to indices of sequences
    (FALSE).}
  \item{\dots}{further arguments to be used by other functions;
    currently not used.}
}
\value{
  The class \code{gengraph} is a list with the following
  components:\cr
  \item{graph}{a graph of class \code{\link[igraph]{igraph}}.}
  \item{clust}{a list containing group information: \code{$membership}:
    an integer giving group membership; \code{$csize}: the size of each
    cluster; \code{$no}: the number of clusters}
  \item{cutoff}{the value used as a cutoff point}
  \item{col}{the color used to plot each group.}
}
\seealso{
  The \code{\link[igraph]{igraph}} package.
}
\author{
  Original idea by Anne Cori and Christophe Fraser.
  Implementation by Thibaut Jombart \email{t.jombart@imperial.ac.uk}.
 }
\examples{
\dontrun{
dat <- haploGen()
res <- gengraph(dat$seq, ngrp=1)
plot(res$graph)
}
}
