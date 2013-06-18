\name{SNPbin-class}
\docType{class}
\alias{SNPbin}
\alias{SNPbin-class}
\alias{[,SNPbin-method}
\alias{[,SNPbin,ANY,ANY-method}
\alias{initialize,SNPbin-method}
\alias{show,SNPbin-method}
\alias{nLoc,SNPbin-method}
\alias{$,SNPbin-method}
\alias{$<-,SNPbin-method}
\alias{names,SNPbin-method}
\alias{ploidy,SNPbin-method}
\alias{ploidy<-,SNPbin-method}
\alias{coerce,SNPbin,integer-method}
\alias{as.integer.SNPbin}
\alias{NA.posi,SNPbin-method}
\alias{cbind.SNPbin}
\alias{c.SNPbin}
\alias{as,integer,SNPbin-method}
\alias{as,numeric,SNPbin-method}


% \alias{,SNPbin-method}
% \alias{,SNPbin-method}
% \alias{,SNPbin-method}
% \alias{,SNPbin-method}
%%%%
\title{Formal class "SNPbin"}
\description{
  The class \code{SNPbin} is a formal (S4) class for storing a genotype
  of binary SNPs in a compact way, using a bit-level coding scheme.
  This storage is most efficient with haploid data, where the memory
  taken to represent data can reduced more than 50 times. However,
  \code{SNPbin} can be used for any level of ploidy, and still remain an
  efficient storage mode.

  A \code{SNPbin} object can be constructed from
  a vector of integers giving the number of the second allele for each
  locus.

  \code{SNPbin} stores a single genotype. To store multiple genotypes,
  use the \linkS4class{genlight} class.
}
\section{Objects from the class SNPbin}{
  \code{SNPbin} objects can be created by calls to \code{new("SNPbin",
    ...)}, where '...' can be the following arguments:
  
  \describe{
    \item{\code{snp}}{a vector of integers or numeric giving numbers of
      copies of the second alleles for each locus. If only one unnamed
      argument is provided to 'new', it is considered as this one.}
    \item{\code{ploidy}}{an integer indicating the ploidy of the
      genotype; if not provided, will be guessed from the data (as the
      maximum from the 'snp' input vector).}
    \item{\code{label}}{an optional character string serving as a label
      for the genotype.}
  }
}
\section{Slots}{
  The following slots are the content of instances of the class
  \code{SNPbin}; note that in most cases, it is better to retrieve
  information via accessors (see below), rather than by accessing the
  slots manually.
  \describe{
    \item{\code{snp}:}{a list of vectors with the class \code{raw}.}
    \item{\code{n.loc}:}{an integer indicating the number of SNPs of the
      genotype.}
    \item{\code{NA.posi}:}{a vector of integer giving the position of
      missing data.}
    \item{\code{label}:}{an optional character string serving as a label
      for the genotype..}
    \item{\code{ploidy}:}{an integer indicating the ploidy of the genotype.}
  }
}
\section{Methods}{
  Here is a list of methods available for \code{SNPbin} objects. Most of
  these methods are accessors, that is, functions which are used to
  retrieve the content of the object. Specific manpages can exist for
  accessors with more than one argument. These are indicated by a '*'
  symbol next to the method's name. This list also contains methods
  for conversion from \code{SNPbin} to other classes.
  \describe{
    \item{[}{\code{signature(x = "SNPbin")}: usual method to subset
      objects in R. The argument indicates how SNPs are to be
      subsetted. It can be a vector of signed integers or of logicals.}
    \item{show}{\code{signature(x = "SNPbin")}: printing of the
      object.}
    \item{$}{\code{signature(x = "SNPbin")}: similar to the @ operator;
      used to access the content of slots of the object.}
    \item{$<-}{\code{signature(x = "SNPbin")}: similar to the @ operator;
      used to replace the content of slots of the object.}
    \item{nLoc}{\code{signature(x = "SNPbin")}: returns the number of
      SNPs in the object.}
    \item{names}{\code{signature(x = "SNPbin")}: returns the names of
      the slots of the object.}
    \item{ploidy}{\code{signature(x = "SNPbin")}: returns the ploidy of
      the genotype.}
    \item{as.integer}{\code{signature(x = "SNPbin")}: converts a
      \code{SNPbin} object to a vector of integers. The S4 method 'as' can
      be used as well (e.g. as(x, "integer")).}
    \item{cbind}{\code{signature(x = "SNPbin")}: merges genotyping of
      the same individual at different SNPs (all stored as
      \linkS4class{SNPbin} objects) into a single \linkS4class{SNPbin}.}
    \item{c}{\code{signature(x = "SNPbin")}: same as cbind.SNPbin.}
  }
}
\author{Thibaut Jombart (\email{t.jombart@imperial.ac.uk})}
\seealso{
  Related class:\cr
  -  \code{\linkS4class{genlight}}, for storing multiple binary SNP
  genotypes. \cr
  -  \code{\linkS4class{genind}}, for storing other types of genetic markers. \cr
}
\examples{
\dontrun{
#### HAPLOID EXAMPLE ####
## create a genotype of 100,000 SNPs
dat <- sample(c(0,1,NA), 1e5, prob=c(.495, .495, .01), replace=TRUE)
dat[1:10]
x <- new("SNPbin", dat)
x
x[1:10] # subsetting
as.integer(x[1:10])

## try a few accessors
ploidy(x)
nLoc(x)
head(x$snp[[1]]) # internal bit-level coding

## check that conversion is OK
identical(as(x, "integer"),as.integer(dat)) # SHOULD BE TRUE

## compare the size of the objects
print(object.size(dat), unit="auto")
print(object.size(x), unit="auto")
object.size(dat)/object.size(x) # EFFICIENCY OF CONVERSION


#### TETRAPLOID EXAMPLE ####
## create a genotype of 100,000 SNPs
dat <- sample(c(0:4,NA), 1e5, prob=c(rep(.995/5,5), 0.005), replace=TRUE)
x <- new("SNPbin", dat)
identical(as(x, "integer"),as.integer(dat)) # MUST BE TRUE

## compare the size of the objects
print(object.size(dat), unit="auto")
print(object.size(x), unit="auto")
object.size(dat)/object.size(x) # EFFICIENCY OF CONVERSION


#### c, cbind ####
a <- new("SNPbin", c(1,1,1,1,1))
b <- new("SNPbin", c(0,0,0,0,0))
a
b
ab <- c(a,b)
ab
identical(c(a,b),cbind(a,b))
as.integer(ab)
}
}
\keyword{classes}
