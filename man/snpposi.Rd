\encoding{UTF-8}
\name{snpposi}

\alias{snpposi.plot}
\alias{snpposi.plot.integer}
\alias{snpposi.plot.numeric}
\alias{snpposi.plot.DNAbin}

\alias{snpposi.test}
\alias{snpposi.test.integer}
\alias{snpposi.test.numeric}
\alias{snpposi.test.DNAbin}

\title{Analyse the position of polymorphic sites}
\description{
  These functions are used to describe the distribution of polymorphic
  sites (SNPs) in an alignment.

  The function \code{snpposi.plot} plots the positions and density of
  SNPs in the alignment.

  The function \code{snpposi.test} tests whether SNPs are randomly
  distributed in the genome, the alternative hypothesis being that they
  are clustered. This test is based on the distances of each SNP to the
  closest SNP. This provides one measure of clustering for each SNP.
  Different statistics can be used to summarise these values (argument
  \code{stat}), but by default the statistics used is the median.

  \code{snpposi.plot} and \code{snpposi.test} are generic functions with
  methods for vectors of integers or numeric (indicating SNP position),
  and for \code{\link[ape]{DNAbin}} objects.
}
\usage{
snpposi.plot(\dots)

\method{snpposi.plot}{integer}(x, genome.size, smooth=0.1,
              col="royalblue", alpha=.2, codon=TRUE, start.at=1, \dots)

\method{snpposi.plot}{numeric}(x, \dots)

\method{snpposi.plot}{DNAbin}(x, \dots)


snpposi.test(\dots)

\method{snpposi.test}{integer}(x, genome.size, n.sim=999, stat=median, \dots)

\method{snpposi.test}{numeric}(x, \dots)

\method{snpposi.test}{DNAbin}(x, \dots)



}
\arguments{
  \item{x}{a vector of integers or numerics containing SNP positions, or
    a set of aligned sequences in a \code{DNAbin} object.}

  \item{genome.size}{an integer indicating the length of genomes.}

  \item{smooth}{a smoothing parameter for the density estimation;
    smaller values will give more local peaks; values have to be
    positive but can be less than 1.}

  \item{col}{the color to be used for the plot; ignored if codon
    positions are represented.}

  \item{alpha}{the alpha level to be used for transparency (density curve).}

  \item{codon}{a logical indicating if codon position should be
    indicated (TRUE, default) or not.}

  \item{start.at}{an integer indicating at which base of a codon the
    alignment starts (defaults to 1); values other than 1, 2 and 3 will
    be ignored.}

  \item{n.sim}{an integer indicating the number of randomizations to be
    used in the Monte Carlo test.}

  \item{stat}{a function used to summarize the measure of physical
    proximity between SNPs; by default, the median is used.}
    \item{\dots}{further arguments to be passed to the \code{integer}
    method.}
}
\value{
  A Monte Carlo test of the class \code{randtest}.
}
\seealso{
  The \code{\link{fasta2DNAbin}} to read fasta alignments with minimum
  RAM use.
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}.
}
\examples{
if(require(ape)){
data(woodmouse)
snpposi.plot(woodmouse, codon=FALSE)
snpposi.plot(woodmouse)

\dontrun{
snpposi.test(c(1,3,4,5), 100)
snpposi.test(woodmouse)
}
}
}

