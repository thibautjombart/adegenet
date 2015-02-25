\encoding{UTF-8}
\name{a-score}
\alias{a.score}
\alias{optim.a.score}
\title{Compute and optimize a-score for Discriminant Analysis of Principal Components (DAPC)}
\description{
  These functions are under development. Please email the author before
  using them for published results.
}
\usage{
a.score(x, n.sim=10, \ldots)

optim.a.score(x, n.pca=1:ncol(x$tab), smart=TRUE, n=10, plot=TRUE,
              n.sim=10, n.da=length(levels(x$grp)), \ldots)
}
\arguments{
\item{x}{a \code{dapc} object.}
\item{n.pca}{a vector of \code{integers} indicating the number of axes retained in the
  Principal Component Analysis (PCA) steps of DAPC. \code{nsim} DAPC
  will be run for each value in \code{n.pca}, unless the smart approach
  is used (see details).}
\item{smart}{a \code{logical} indicating whether a smart, less
  computer-intensive approach should be used (TRUE, default) or not
  (FALSE). See details section.}
\item{n}{an \code{integer} indicating the numbers of values spanning the
  range of \code{n.pca} to be used in the smart approach.}
\item{plot}{a \code{logical} indicating whether the results should be
  displayed graphically (TRUE, default) or not (FALSE).}
\item{n.sim}{an \code{integer} indicating the number of simulations to
  be performed for each number of retained PC.}
\item{n.da}{an \code{integer} indicating the number of axes retained in the
  Discriminant Analysis step.}
\item{\ldots}{further arguments passed to other methods; currently unused..}
}
\details{
  The Discriminant Analysis of Principal Components seeks a reduced
  space inside which observations are best discriminated into
  pre-defined groups. One way to assess the quality of the
  discrimination is looking at re-assignment of individuals to their
  prior group, successful re-assignment being a sign of strong discrimination.

  However, when the original space is very large, ad hoc solutions can
  be found, which discriminate very well the sampled individuals but
  would perform poorly on new samples. In such a case, DAPC
  re-assignment would be high even for randomly chosen clusters.
  The a-score measures this bias. It is computed as (Pt-Pr), where Pt is
  the reassignment probability using the true cluster, and Pr is the
  reassignment probability for randomly permuted clusters. A a-score
  close to one is a sign that the DAPC solution is both strongly
  discriminating and stable, while low values (toward 0 or lower)
  indicate either weak discrimination or instability of the results.

  The a-score can serve as a criterion for choosing the optimal number of
  PCs in the PCA step of DAPC, i.e. the number of PC maximizing the
  a-score. Two procedures are implemented in \code{optim.a.score}. The
  smart procedure selects evenly distributed number of PCs in a
  pre-defined range, compute the a-score for each, and then interpolate
  the results using splines, predicting an approximate optimal number of
  PCs. The other procedure (when \code{smart} is FALSE) performs the
  computations for all number of PCs request by the user. The 'optimal'
  number is then the one giving the highest mean a-score (computed over
  the groups).
}
\value{
  === a.score ===\cr
 \code{a.score} returns a list with the following components:\cr
  \item{tab}{a matrix of a-scores with groups in columns and simulations in row.}
  \item{pop.score}{a vector giving the mean a-score for each population.}
  \item{mean}{the overall mean a-score.}\cr
  
  === optim.a.score ===\cr
  \code{optima.score} returns a list with the following components:\cr
   \item{pop.score}{a list giving the mean a-score of the populations
   for each number of retained PC (each element of the list corresponds
   to a number of retained PCs).}
  \item{mean}{a vector giving the overall mean a-score for each number
    of retained PCs.}
  \item{pred}{(only when \code{smart} is TRUE) the predictions of the
    spline, given in x and y coordinates.}
 \item{best}{the optimal number of PCs to be retained.}
}
\references{
  Jombart T, Devillard S and Balloux F  (2010) Discriminant analysis of
  principal components: a new method for the analysis of genetically
  structured populations. BMC Genetics11:94. doi:10.1186/1471-2156-11-94
}
\seealso{
  - \code{\link{find.clusters}}: to identify clusters without prior.

  - \code{\link{dapc}}: the Discriminant Analysis of Principal
    Components (DAPC)
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\keyword{multivariate}