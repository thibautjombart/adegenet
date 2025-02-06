\encoding{UTF-8}
\name{pairDistPlot}
\alias{pairDist}
\alias{pairDist.default}
\alias{pairDistPlot}
\alias{pairDistPlot.default}
\alias{pairDistPlot.dist}
\alias{pairDistPlot.matrix}
\alias{pairDistPlot.genind}
\alias{pairDistPlot.DNAbin}
\title{Pairwise distance plots}
\description{
  The function \code{pairDistPlot} extracts and plots pairwise distances
  between different groups (graphs use ggplot2). The function
  \code{pairDistPlot} does the same, without the graphs.
  
  \code{pairDistPlot} is a generic function with methods for the
 following types of objects:\cr
 - \code{dist} \cr
 - \code{matrix} (only numeric data)\cr
 - \code{\linkS4class{genind}} objects (genetic markers, individuals)\cr
 -  \code{\link[ape]{DNAbin}} objects (DNA sequences)
}
\usage{
pairDist(x, \dots)

pairDistPlot(x, \dots)

\method{pairDistPlot}{dist}(x, grp, within=FALSE, sep="-", data=TRUE,
             violin=TRUE, boxplot=TRUE, jitter=TRUE, \dots)

\method{pairDistPlot}{matrix}(x, grp, within=FALSE, sep="-", data=TRUE,
             violin=TRUE, boxplot=TRUE, jitter=TRUE, \dots)

\method{pairDistPlot}{genind}(x, grp, within=FALSE, sep="-", data=TRUE,
             violin=TRUE, boxplot=TRUE, jitter=TRUE, \dots)

\method{pairDistPlot}{DNAbin}(x, grp, within=FALSE, sep="-", data=TRUE,
             violin=TRUE, boxplot=TRUE, jitter=TRUE, \dots)

}
\arguments{
  \item{x}{pairwise distances provided as a \code{dist} or a symmetric
    \code{matrix}, or \code{\linkS4class{genind}} or \code{DNAbin}
    object. For \code{\linkS4class{genind}} objects, pairwise squared
    Euclidean distances are computed from the allele data. For
    \code{DNAbin} objects, distances are computed uing
    \code{dist.dna}, and '...' is used to pass arguments to the
    function.  }

  \item{grp}{a factor defining a grouping of individuals.}

  \item{within}{a logical indicating whether to keep within-group
    comparisons.}
  
  \item{sep}{a character used as separator between group names}

  \item{data}{a logical indicating whether data of the plot should be returned.}

  \item{violin}{a logical indicating whether a violinplot should be generated.}

  \item{boxplot}{a logical indicating whether a boxplot should be generated.}

  \item{jitter}{a logical indicating whether a jitter-plot should be generated.}

  \item{\dots}{further arguments to be used by other functions; used for
    \code{DNAbin} object to pass argumetns to \code{dist.dna}.}
}
\value{
  A list with different components, depending on the values of the
    arguments. Plots are returned as \code{ggplot2} objects.
}
\seealso{
  \code{\link{gengraph}} to identify connectivity based on distances.
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}.
 }
\examples{
\dontrun{

## use a subset of influenza data
data(H3N2)
set.seed(1)
dat <- H3N2[sample(1:nInd(H3N2), 100)]

## get pairwise distances
temp <- pairDistPlot(dat, other(dat)$epid)

## see raw data
head(temp$data)

## see plots
temp$boxplot
temp$violin
temp$jitter

}
}
