\encoding{UTF-8}
\name{DAPC cross-validation}
\alias{xvalDapc}
\alias{xvalDapc.default}
\alias{xvalDapc.data.frame}
\alias{xvalDapc.matrix}
\alias{xvalDapc.genlight}
\alias{xvalDapc.genind}
\title{Cross-validation for Discriminant Analysis of Principal Components (DAPC)}
\description{
  The function \code{xvalDapc} performs stratified cross-validation of DAPC 
  using varying numbers of PCs (and keeping the number of discriminant
  functions fixed); \code{xvalDapc} is a generic with methods for
  \code{data.frame} and \code{matrix}.\cr
}
\usage{
xvalDapc(x, \dots)

\method{xvalDapc}{default}(x, grp, n.pca.max = 300, n.da = NULL,
              training.set = 0.9, result = c("groupMean", "overall"),
              center = TRUE, scale = FALSE,
              n.pca=NULL, n.rep = 30, xval.plot = TRUE, \dots)

\method{xvalDapc}{data.frame}(x, grp, n.pca.max = 300, n.da = NULL,
              training.set = 0.9, result = c("groupMean", "overall"),
              center = TRUE, scale = FALSE,
              n.pca=NULL, n.rep = 30, xval.plot = TRUE, \dots)

\method{xvalDapc}{matrix}(x, grp, n.pca.max = 300, n.da = NULL,
              training.set = 0.9, result = c("groupMean", "overall"),
              center = TRUE, scale = FALSE,
              n.pca=NULL, n.rep = 30, xval.plot = TRUE, \dots)

\method{xvalDapc}{genlight}(x, \dots)

\method{xvalDapc}{genind}(x, \dots)
}
\arguments{
  \item{x}{A numeric \code{data.frame} or \code{matrix}, or a
    \linkS4class{genind} or \linkS4class{genlight} object. Values must be
    finite and complete; missing values are not imputed. The genetic-object
    methods convert to allele counts or SNP dosages, respectively.}
  \item{grp}{a \code{factor} indicating the group membership of
    individuals.}
  \item{n.pca.max}{maximum number of PCA components to retain.}
  \item{n.da}{an \code{integer} indicating the number of axes retained in the
    Discriminant Analysis step. If \code{NULL}, n.da defaults to 1 less than
    the number of groups.}
  \item{training.set}{the proportion of data (individuals) to be used
    for the training set; defaults to 0.9 if all groups have >= 10 members; otherwise,
    training.set scales automatically to the largest proportion that still 
    ensures all groups will be present in both training and validation sets.}
  \item{result}{a character string; "groupMean" for group-wise assignment
    sucess, or  "overall" for an overall mean assignment success; see details.}
  \item{center}{a \code{logical} indicating whether variables should be centred to
    mean 0 (TRUE, default) or not (FALSE). Centring is fitted separately
    within each training split.}
  \item{scale}{a \code{logical} indicating whether variables should be scaled
    (TRUE) or not (FALSE, default). Scaling consists in dividing variables by their
    (estimated) standard deviation to account for trivial differences in
    variances.}
  \item{n.pca}{an \code{integer} vector indicating the number of
    different number of PCA axes to be retained for the cross
    validation; if \code{NULL}, this will be dertermined automatically.}
  \item{n.rep}{the number of replicates to be carried out at each level of
    PC retention; defaults to 30.}
  \item{xval.plot}{a logical indicating whether a plot of the cross-validation
    results should be generated.}
  \item{\dots}{further arguments to be passed to \code{\link[boot]{boot}}.
    \strong{see Details.}} }
\details{ The Discriminant Analysis of Principal Components (DAPC) relies on
  dimension reduction of the data using PCA followed by a linear discriminant
  analysis. How many PCA axes to retain is often a non-trivial question. Cross
  validation provides an objective way to decide how many axes to retain:
  different numbers are tried and the quality of the corresponding DAPC is
  assessed by cross- validation: DAPC is performed on a training set, typically
  made of 90\% of the observations (comprising 90\% of the observations in each
  subpopulation) , and then used to predict the groups of the 10\% of remaining
  observations. The current method uses the average prediction success per group
  (result="groupMean"), or the overall prediction success (result="overall").
  The number of PCs associated with the lowest Mean Squared Error is then
  retained in the DAPC.
  \subsection{Training-only preprocessing}{
    Each replicate fits PCA, including centring and scaling, using only
    the training individuals. Validation individuals are projected using
    that training transformation before their memberships are predicted.
    No PCA fitted to the full dataset is used to score validation samples.
    The final returned DAPC is then fitted to all input individuals, using
    the selected PC count and the requested \code{center}, \code{scale},
    and \code{n.da} settings.

    A PC count unsupported by the PCA rank or the pooled within-group
    LDA rank of any training split is excluded from selection, with a
    warning. The latter is at most \eqn{N-K}, for \eqn{N} training
    individuals in \eqn{K} represented groups. Such replicates have
    \code{NA} success, and the candidate's mean success and RMSE are
    \code{NA}; they are not silently evaluated using fewer PCs.
    If no candidate remains, reduce the requested PC counts.
    Singleton groups remain in training and cannot be evaluated in
    validation. Every other group contributes at least one individual to
    each set.
  }
  \subsection{Missing data and interpretation}{
    Missing values cause an explicit error. Unlike \code{dapc.genind},
    this function does not automatically replace them with allele means.
    Imputing the full dataset before cross-validation can allow validation
    individuals to influence training preprocessing. To evaluate an
    imputation-plus-DAPC workflow, imputation must also be fitted within
    each training split; this function does not currently implement that
    workflow. Mean replacement does not propagate genotype uncertainty.

    These cross-validation results select the PC count; they are not an
    independent performance estimate for that selected model. Independent
    test data or an outer validation loop are needed for that purpose.
    Splits are stratified by group, not blocked by sequencing batch,
    family, or collection site. Strong assignment alone therefore does
    not rule out technical or sampling confounding.
  }
  \subsection{Parallel Computing}{
    The permutation of the data for cross-validation is performed in part by the
    function\code{\link[boot]{boot}}. If you have a modern computer, it is
    likely that you have multiple cores on your system. R by default utilizes
    only one of these cores unless you tell it otherwise. For details, please
    see the documentation of \code{\link[boot]{boot}}. Basically, if you want to
    use multiple cores, you need two arguments:
    \enumerate{
      \item \code{parallel} - what R parallel system to use (see below)
      \item \code{ncpus} - number of cores you want to use
    } 
    If you are on a unix system (Linux or OSX), you will want to specify
    \code{parallel = "multicore"}. If you are on Windows, you will want to
    specify \code{parallel = "snow"}.
  } 
}
\value{
 A \code{list} containing seven items, and a \code{plot} of the results. The
 first is a \code{data.frame} with two columns, the first giving the number of
 PCs of PCA retained in the corresponding DAPC, and the second giving the
 proportion of successful group assignment for each replicate. The second item
 gives the mean and confidence interval for random chance. The third gives the
 mean successful assignment at each level of PC retention. The fourth indicates
 which number of PCs is associated with the highest mean success. The fifth
 gives the Root Mean Squared Error at each level of PC retention. The sixth
 indicates which number of PCs is associated with the lowest MSE. The seventh
 item contains the DAPC carried out with the optimal number of PCs, determined
 with reference to MSE.
 
 If \code{xval.plot=TRUE} a scatterplot of the results of cross-validation 
 will be displayed. 
 
}
\references{
  Jombart T, Devillard S and Balloux F  (2010) Discriminant analysis of
  principal components: a new method for the analysis of genetically
  structured populations. BMC Genetics11:94. doi:10.1186/1471-2156-11-94
}
\seealso{
  \code{\link{dapc}}
}
\author{ Caitlin Collins \email{caitlin.collins12@imperial.ac.uk},
Thibaut Jombart \email{t.jombart@imperial.ac.uk},
Zhian N. Kamvar \email{kamvarz@science.oregonstate.edu} }
\examples{
\dontrun{
## CROSS-VALIDATION ##
data(sim2pop)
xval <- xvalDapc(sim2pop@tab, pop(sim2pop), n.pca.max=100, n.rep=3)
xval

## 100 replicates ##

# Serial version (SLOW!)
system.time(xval <- xvalDapc(sim2pop@tab, pop(sim2pop), n.pca.max=100, n.rep=100))

# Parallel version (faster!)
system.time(xval <- xvalDapc(sim2pop@tab, pop(sim2pop), n.pca.max=100, n.rep=100, 
                             parallel = "multicore", ncpus = 2))
}
}
\keyword{multivariate}
