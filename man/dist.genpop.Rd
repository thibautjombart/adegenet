\encoding{UTF-8}
\name{dist.genpop}
\alias{dist.genpop}
\title{ Genetic distances between populations }
\description{
  This function computes measures of genetic distances
  between populations using a \code{genpop} object. \cr
  Currently, five distances are available, some of which are euclidian
  (see details).\cr
  
  A non-euclidian distance can be transformed into an Euclidean one
  using \code{\link[ade4]{cailliez}} in order to perform a
  Principal Coordinate Analysis \code{\link[ade4]{dudi.pco}} (both
  functions in \code{ade4}). \cr
  
  The function \code{dist.genpop} is based on former \code{dist.genet}
  function of \code{ade4} package.
}
\usage{
dist.genpop(x, method = 1, diag = FALSE, upper = FALSE)
}
\arguments{
  \item{x}{ a list of class \code{genpop}  }
  \item{method}{ an integer between 1 and 5. See details }
  \item{diag}{ a logical value indicating whether the diagonal of the distance matrix should be printed by \code{print.dist} }
  \item{upper}{ a logical value indicating whether the upper triangle of the distance matrix should be printed by \code{print.dist} }
}
\details{
Let \bold{A} a table containing allelic frequencies with \emph{t} populations (rows) and \emph{m} alleles (columns).\cr
Let  \eqn{\nu} the number of loci. The locus \emph{j} gets \emph{m(j)} alleles. 
\eqn{m=\sum_{j=1}^{\nu} m(j)}\cr

For the row \emph{i} and the modality \emph{k} of the variable \emph{j}, notice the value \eqn{a_{ij}^k} (\eqn{1 \leq i \leq t}, \eqn{1 \leq j \leq \nu},
\eqn{1 \leq k \leq m(j)}) the value of the initial table.\cr

\eqn{a_{ij}^+=\sum_{k=1}^{m(j)}a_{ij}^k} and \eqn{p_{ij}^k=\frac{a_{ij}^k}{a_{ij}^+}}\cr

Let \bold{P} the table of general term \eqn{p_{ij}^k}\cr
\eqn{p_{ij}^+=\sum_{k=1}^{m(j)}p_{ij}^k=1}, \eqn{p_{i+}^+=\sum_{j=1}^{\nu}p_{ij}^+=\nu}, \eqn{p_{++}^+=\sum_{j=1}^{\nu}p_{i+}^+=t\nu}\cr

The option \code{method} computes the distance matrices between populations using the frequencies \eqn{p_{ij}^k}. \cr

1. Nei's distance (not Euclidean): \cr
\eqn{D_1(a,b)=- \ln(\frac{\sum_{k=1}^{\nu} \sum_{j=1}^{m(k)}
p_{aj}^k p_{bj}^k}{\sqrt{\sum_{k=1}^{\nu} \sum_{j=1}^{m(k)}
{(p_{aj}^k) }^2}\sqrt{\sum_{k=1}^{\nu} \sum_{j=1}^{m(k)}
{(p_{bj}^k)}^2}})}\cr

2. Angular distance or Edwards' distance (Euclidean):\cr
\eqn{D_2(a,b)=\sqrt{1-\frac{1}{\nu} \sum_{k=1}^{\nu}
\sum_{j=1}^{m(k)} \sqrt{p_{aj}^k  p_{bj}^k}}}\cr

3. Coancestrality coefficient or Reynolds' distance (Eucledian):\cr
\eqn{D_3(a,b)=\sqrt{\frac{\sum_{k=1}^{\nu}
\sum_{j=1}^{m(k)}{(p_{aj}^k - p_{bj}^k)}^2}{2 \sum_{k=1}^{\nu} (1-
\sum_{j=1}^{m(k)}p_{aj}^k p_{bj}^k)}}}\cr

4. Classical Euclidean distance or Rogers' distance (Eucledian):\cr
\eqn{D_4(a,b)=\frac{1}{\nu} \sum_{k=1}^{\nu} \sqrt{\frac{1}{2}
\sum_{j=1}^{m(k)}{(p_{aj}^k - p_{bj}^k)}^2}}\cr

5. Absolute genetics distance or Provesti 's distance (not Euclidean):\cr
\eqn{D_5(a,b)=\frac{1}{2{\nu}} \sum_{k=1}^{\nu} \sum_{j=1}^{m(k)}
|p_{aj}^k - p_{bj}^k|}
}
\value{
returns a distance matrix of class \code{dist} between the rows of the data frame
}
\references{ 
To complete informations about distances:\cr

Distance 1:\cr
Nei, M. (1972) Genetic distances between populations. \emph{American Naturalist}, \bold{106}, 283--292. \cr
Nei M. (1978) Estimation of average heterozygosity and genetic distance from a small number of individuals. \emph{Genetics}, \bold{23}, 341--369. \cr
Avise, J. C. (1994) Molecular markers, natural history and evolution. Chapman & Hall, London. 

Distance 2:\cr
Edwards, A.W.F. (1971) Distance between populations on the basis of gene frequencies. \emph{Biometrics}, \bold{27}, 873--881. \cr
Cavalli-Sforza L.L. and Edwards A.W.F. (1967) Phylogenetic analysis: models and estimation procedures. \emph{Evolution}, \bold{32}, 550--570. \cr
Hartl, D.L. and Clark, A.G. (1989) Principles of population genetics. Sinauer Associates, Sunderland, Massachussetts (p. 303).

Distance 3:\cr
Reynolds, J. B., B. S. Weir, and C. C. Cockerham. (1983) Estimation of the coancestry coefficient: basis for a short-term genetic distance. \emph{Genetics}, \bold{105}, 767--779.

Distance 4:\cr
Rogers, J.S. (1972) Measures of genetic similarity and genetic distances. \emph{Studies in Genetics}, Univ. Texas Publ., \bold{7213}, 145--153.  \cr
Avise, J. C. (1994) Molecular markers, natural history and evolution. Chapman & Hall, London. 

Distance 5:\cr
Prevosti A. (1974) La distancia genetica entre poblaciones. \emph{Miscellanea Alcobe}, \bold{68}, 109--118. \cr
Prevosti A., Oca\~na J. and Alonso G. (1975) Distances between populations of Drosophila subobscura, based on chromosome arrangements frequencies. \emph{Theoretical and Applied Genetics}, \bold{45}, 231--241. \cr

For more information on dissimilarity indexes:\cr
Gower J. and Legendre P. (1986) Metric and Euclidean properties of
dissimilarity coefficients. \emph{Journal of Classification}, \bold{3},
5--48 \cr

Legendre P. and Legendre L. (1998) \emph{Numerical Ecology}, Elsevier
Science B.V. 20, pp274--288.\cr 
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}\cr
Former dist.genet code by Daniel Chessel \email{chessel@biomserv.univ-lyon1.fr}\cr
and documentation by Anne B. Dufour \email{dufour@biomserv.univ-lyon1.fr}
}
\seealso{ \code{\link[ade4]{cailliez}},\code{\link[ade4]{dudi.pco}} 
}
\examples{
\dontrun{
data(microsatt)
obj <- as.genpop(microsatt$tab)

listDist <- lapply(1:5, function(i) cailliez(dist.genpop(obj,met=i)))
for(i in 1:5) {attr(listDist[[i]],"Labels") <- obj@pop.names}
listPco <- lapply(listDist, dudi.pco,scannf=FALSE)

par(mfrow=c(2,3))
for(i in 1:5) {scatter(listPco[[i]],sub=paste("Dist:", i))}

}
}
\keyword{multivariate}
