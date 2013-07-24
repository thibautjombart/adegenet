\encoding{UTF-8}
\name{DAPC cross-validation}
\alias{xvalDapc}
\alias{xvalDapc.data.frame}
\alias{xvalDapc.matrix}
\title{Cross-validation for Discriminant Analysis of Principal Components (DAPC)}
\description{
  The function \code{xvalDapc} performs stratified cross-validation of DAPC 
  using varying numbers of PCs (and keeping the number of discriminant
  functions fixed); \code{xvalDapc} is a generic with methods for
  \code{data.frame} and \code{matrix}.\cr
}
\usage{
xvalDapc(x, \dots)

\method{xvalDapc}{data.frame}(x, grp, n.pca.max, n.da=NULL,
              training.set = 0.9, result=c("groupMean","overall"),
              center=TRUE, scale=FALSE,
              n.pca=NULL, n.rep=10, \dots)

\method{xvalDapc}{matrix}(x, grp, n.pca.max, n.da=NULL,
              training.set = 0.9, result=c("groupMean","overall"),
              center=TRUE, scale=FALSE,
              n.pca=NULL, n.rep=10, \dots)
}
\arguments{
  \item{x}{\code{a data.frame} or a \code{matrix} used as input of DAPC.}
  \item{grp}{a \code{factor} indicating the group membership of
    individuals.}
  \item{n.pca.max}{maximum number of PCA components to retain.}
  \item{n.da}{an \code{integer} indicating the number of axes retained in the
    Discriminant Analysis step. If \code{NULL}, interactive selection is triggered.}
  \item{training.set}{the proportion of data (individuals) to be used
    for the training set; defaults to 0.9, if all groups have >= 10 members; otherwise,
    training.set scales automatically to the largest proportion that still 
    ensures all groups will be present in both training and validation sets.}
  \item{result}{a character string; "groupMean" for group-wise assignment
    sucess, or  "overall" for an overall mean assignment success; see details.}
  \item{center}{a \code{logical} indicating whether variables should be centred to
    mean 0 (TRUE, default) or not (FALSE). Always TRUE for \linkS4class{genind} objects.}
  \item{scale}{a \code{logical} indicating whether variables should be scaled
    (TRUE) or not (FALSE, default). Scaling consists in dividing variables by their
    (estimated) standard deviation to account for trivial differences in
    variances.}
  \item{n.pca}{an \code{integer} vector indicating the number of
    different number of PCA axes to be retained for the cross
    validation; if \code{NULL}, this will be dertermined automatically.}
  \item{n.rep}{the number of replicate to be used for each number of PCA
    components retained.}
  \item{\dots}{further arguments to be passed to other methods.}
}
\details{
  The Discriminant Analysis of Principal Components (DAPC) relies on
  dimension reduction of the data using PCA followed by a linear
  discriminant analysis. How many PCA axes to retain is often a
  non-trivial question. Cross validation provides an objective way to
  decide how many axes to retain: different numbers are tried and the
  quality of the corresponding DAPC is assessed by cross-validation:
  DAPC is performed on a training set, typically made of 90\% of the
  observations (comprising 90\% of the observations in each subpopulation)
  , and then used to predict the groups of the 10\% of remaining observations. 
  The current method uses the average prediction success per group (result="groupMean"), 
  or the overall prediction success (result="overall").
}
\value{
 A \code{data.frame} with two columns, the first giving the proportion
 of successful group assignment, and the second giving the number of PCA
 retained for the corresponding DAPC.
}
\references{
  Jombart T, Devillard S and Balloux F  (2010) Discriminant analysis of
  principal components: a new method for the analysis of genetically
  structured populations. BMC Genetics11:94. doi:10.1186/1471-2156-11-94
}
\seealso{
  \code{\link{dapc}}
}
\authors{ Thibaut Jombart \email{t.jombart@imperial.ac.uk}, 
Caitlin Collins \email{caitlin.collins12@imperial.ac.uk} }
\examples{
\dontrun{
## CROSS-VALIDATION ##
data(sim2pop)
xval <- xvalDapc(sim2pop@tab, pop(sim2pop), n.pca.max=100, n.rep=3)
xval
boxplot(xval$success~xval$n.pca, xlab="Number of PCA components",
ylab="Classification succes", main="DAPC - cross-validation")

}
}
\keyword{multivariate}
