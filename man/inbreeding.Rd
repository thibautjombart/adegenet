\encoding{UTF-8}
\name{Inbreeding estimation}
\alias{inbreeding}
\title{Likelihood-based estimation of inbreeding}
\description{
  The function \code{inbreeding} estimates the inbreeding coefficient
  of an individuals (F) by computing its likelihood function. It can
  return either the density of probability of F, or a sample of F values
  from this distribution. This operation is performed for all the
  individuals of a \linkS4class{genind} object. Any ploidy greater than
  1 is acceptable.
}
\usage{
inbreeding(x, pop = NULL, truenames = TRUE, res.type = c("sample","function"),
           N = 200, M = N * 10)
}
\arguments{
  \item{x}{an object of class \linkS4class{genind}.}
  \item{pop}{a factor giving the 'population' of each individual. If NULL,
    pop is seeked from \code{pop(x)}. Note that the term population refers in
    fact to any grouping of individuals'.}
  \item{truenames}{a logical indicating whether true names should be
    used (TRUE, default) instead of generic labels (FALSE); used if
    res.type is "matrix".}
 \item{res.type}{a character string matching "sample" or "function",
   specifying whether the output should be a function giving the density
   of probability of F values ("function") or a sample of F values taken
   from this distribution ("sample", default).}
 \item{N}{an integer indicating the size of the sample to be taken from
   the distribution of F values.}
 \item{M}{an integer indicating the number of different F values to be
   used to generate the sample. Values larger than N are recommended to
   avoid poor sampling of the distribution.}
}
\value{
  A named list with one component for each individual, each of which is
  a function or a vector of sampled F values (see \code{res.type} argument).
}
\seealso{
  \code{\link{Hs}}: computation of expected heterozygosity.
}
\details{
  Let \eqn{F} denote the inbreeding coefficient, defined as the
  probability for an individual to inherit two identical alleles from a
  single ancestor.

  Let \eqn{p_i} refer to the frequency of allele \eqn{i} in the
  population. Let \eqn{h} be an variable which equates 1 if the
  individual is homozygote, and 0 otherwise. For one locus, the
  probability of being homozygote is computed as:

  \eqn{ F + (1-F) \sum_i p_i^2}

  The probability of being heterozygote is:
  \eqn{1 - (F + (1-F) \sum_i p_i^2)}

  The likelihood of a genotype is defined as the probability of being
  the observed state (homozygote or heterozygote). In the case of
  multilocus genotypes, log-likelihood are summed over the loci.
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}\cr
}
\examples{
\dontrun{
## cattle breed microsatellite data
data(microbov)

## isolate Lagunaire breed
lagun <- seppop(microbov)$Lagunaire

## estimate inbreeding - return sample of F values
Fsamp <- inbreeding(lagun, N=30)

## plot the first 10 results
invisible(sapply(Fsamp[1:10], function(e) plot(density(e), xlab="F",
xlim=c(0,1), main="Density of the sampled F values")))

## compute means for all individuals
Fmean=sapply(Fsamp, mean)
hist(Fmean, col="orange", xlab="mean value of F",
main="Distribution of mean F across individuals")

## estimate inbreeding - return proba density functions
Fdens <- inbreeding(lagun, res.type="function")

## view function for the first individual
Fdens[[1]]

## plot the first 10 functions
invisible(sapply(Fdens[1:10], plot, ylab="Density",
main="Density of probability of F values"))
}
}