\encoding{UTF-8}
\name{fasta2genlight}
\alias{fasta2genlight}
\title{ Extract Single Nucleotide Polymorphism (SNPs) from alignments}
\description{
  The function \code{fasta2genlight} reads alignments with the fasta
  format (extensions ".fasta", ".fas", or ".fa"), extracts the binary
  SNPs, and converts the output into a \linkS4class{genlight} object. 

  The function reads data by chunks of a few genomes (minimum 1, no
  maximum) at a time, which allows one to read massive datasets with
  negligible RAM requirements (albeit at a cost of computational
  time). The argument \code{chunkSize} indicates the number of genomes
  read at a time. Increasing this value decreases the computational time
  required to read data in, while increasing memory requirements.

  Multiple cores can be used to decrease the overall computational time
  on parallel architectures (needs the package \code{parallel}).
}
\usage{
fasta2genlight(file, quiet=FALSE, chunkSize = 1000, saveNbAlleles=FALSE,
               parallel = require("parallel"), n.cores = NULL, \dots)
}
\arguments{
  \item{file}{ a character string giving the path to the file to
    convert, with the extension ".fa", ".fas", or ".fasta".}
  \item{quiet}{ logical stating whether a conversion messages should be
    printed (FALSE,default) or not (TRUE).}
  \item{chunkSize}{an integer indicating the number of genomes to be
    read at a time; larger values require more RAM but decrease the time
    needed to read the data.}
  \item{saveNbAlleles}{a logical indicating whether the number of
    alleles for each loci in the original alignment should be saved in the
    \code{other} slot (TRUE), or not (FALSE, default). In large genomes,
    this takes some space but allows for tracking SNPs with more than 2
    alleles, lost during the conversion.}
  \item{parallel}{a logical indicating whether multiple cores -if
    available- should be used for the computations (TRUE, default), or
    not (FALSE); requires the package \code{parallel} to be installed
    (see details).}
  \item{n.cores}{if \code{parallel} is TRUE, the number of cores to be
    used in the computations; if NULL, then the maximum number of cores
    available on the computer is used.}
  \item{\dots}{other arguments to be passed to other functions -
    currently not used.}
}
\details{

  === Using multiple cores ===

  Most recent machines have one or several processors with multiple
  cores. R processes usually use one single core. The package
  \code{parallel} allows for parallelizing some computations on
  multiple cores, which decreases drastically computational time.

  To use this functionality, you need to have the last version of the
  \code{parallel} package installed.
}
\value{an object of the class \linkS4class{genlight}}
\seealso{
  - \code{?genlight} for a description of the class \linkS4class{genlight}.

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
obj <- fasta2genlight(myPath, chunk=10) # process 10 sequences at a time
obj

## look at extracted information
position(obj)
alleles(obj)
locNames(obj)

## plot positions of polymorphic sites
temp <- density(position(obj), bw=10)
plot(temp, xlab="Position in the alignment", lwd=2, main="Location of the SNPs")
points(position(obj), rep(0, nLoc(obj)), pch="|", col="red")
}
}
\keyword{manip}
