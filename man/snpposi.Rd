\encoding{UTF-8}
\name{snpposi}
\alias{snpposi}
\alias{snpposi.DNAbin}
\title{Analyse the position of polymorphic sites}
\description{
  The function \code{snpposi.test} tests whether SNPs are randomly
  distributed in the genome, the alternative hypothesis being that they
  are clustered. This test is based on the distances of each SNP to the
  closest SNP. This provides one measure of clustering for each SNP.
  Different statistics can be used to summarise these values (argument
  \code{stat}), but by default the statistics used is the
  median. \code{snpposi.test} is a generic with methods for vectors of
  integers or numeric (indicating SNP position), and for
  \code{\link[ape]{DNAbin}} objects.
}
\usage{
snpposi(\dots)

\method{snpposi}{integer}(x, genome.size, n.sim=999, stat=median\dots)

\method{snpposi}{numeric}(x, \dots)

\method{snpposi}{DNAbin}(x, \dots)



}
\arguments{
  \item{x}{a vector of integers or numerics containing SNP positions, or
    a set of aligned sequences in a \code{DNAbin} object.}

  \item{genome.size}{an integer indicating the length of genomes.}

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
\dontrun{
data(woodmouse)

snpposi.test(woodmouse)

snpposi.test(c(1,3,4,5), 100)
}
}
