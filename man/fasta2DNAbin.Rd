\encoding{UTF-8}
\name{fasta2DNAbin}
\alias{fasta2DNAbin}
\title{ Read large DNA alignments into R}
\description{
  The function \code{fasta2DNAbin} reads alignments with the fasta
  format (extensions ".fasta", ".fas", or ".fa"), and outputs a
  \code{\link[ape]{DNAbin}} object (the efficient DNA representation from the
  ape package). The output contains either the full alignments, or only
  SNPs. This implementation is designed for memory-efficiency,
  and can read in larger datasets than Ape's \code{\link[ape]{read.dna}}.

  The function reads data by chunks of a few genomes (minimum 1, no
  maximum) at a time, which allows one to read massive datasets with
  negligible RAM requirements (albeit at a cost of computational
  time). The argument \code{chunkSize} indicates the number of genomes
  read at a time. Increasing this value decreases the computational time
  required to read data in, while increasing memory requirements.

}
\usage{
fasta2DNAbin(file, quiet=FALSE, chunkSize=10, snpOnly=FALSE)
}
\arguments{
  \item{file}{ a character string giving the path to the file to
    convert, with the extension ".fa", ".fas", or ".fasta".}
  \item{quiet}{a logical stating whether a conversion messages should be
    printed (FALSE, default) or not (TRUE).}
  \item{chunkSize}{an integer indicating the number of genomes to be
    read at a time; larger values require more RAM but decrease the time
    needed to read the data.}
  \item{snpOnly}{a logical indicating whether SNPs only should be returned.}
}
\value{an object of the class \code{\link[ape]{DNAbin}}}
\seealso{
  - \code{?DNAbin} for a description of the class \code{\link[ape]{DNAbin}}.

  - \code{\link{read.snp}}: read SNPs in adegenet's '.snp' format.

  - \code{\link{read.PLINK}}: read SNPs in PLINK's '.raw' format.

  - \code{\link{df2genind}}: convert any multiallelic markers into
  adegenet \linkS4class{genind}.

  - \code{\link{import2genind}}: read multiallelic markers from various
  software into adegenet.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
## show the example file ##
## this is the path to the file:
myPath <- system.file("files/usflu.fasta",package="adegenet")
myPath

## read the file
obj <- fasta2DNAbin(myPath, chunk=10) # process 10 sequences at a time
obj
}
}
\keyword{manip}
