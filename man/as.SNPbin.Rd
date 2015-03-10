\name{as.SNPbin}
\alias{as,SNPbin,integer-method}
\alias{as,SNPbin,numeric-method}
\alias{as.SNPbin}
\alias{as.SNPbin,integer-method}
\alias{as.SNPbin,numeric-method}
\alias{coerce,integer,SNPbin-method}
\alias{coerce,numeric,SNPbin-method}

% \alias{,SNPbin-method}
% \alias{,SNPbin-method}
% \alias{,SNPbin-method}
% \alias{,SNPbin-method}
%%%%
\title{Conversion to class "SNPbin"}
\description{
  The class \linkS4class{SNPbin} is a formal (S4) class for storing a genotype
  of binary SNPs in a compact way, using a bit-level coding scheme. New
  instances of this class are best created using \code{new}; see the
  manpage of \linkS4class{SNPbin} for more information on this point.

  As a shortcut, conversion methods can be used to convert various
  objects into a \linkS4class{SNPbin} object. Conversions can be
  achieved using S3-style (\code{as.SNPbin(x)}) or S4-style
  (\code{as(x,"SNPbin"}) procedures. All of them call upon the
  constructor (\code{new}) of \linkS4class{SNPbin} objects.

  Conversion is currently available from the following objects:
  - integer vectors
  - numeric vectors
}
\author{Thibaut Jombart (\email{t.jombart@imperial.ac.uk})}
\seealso{
  Related class:\cr
  -  \code{\linkS4class{SNPbin}}
  -  \code{\linkS4class{genlight}}, for storing multiple binary SNP
  genotypes. \cr
}
\examples{
\dontrun{
## data to be converted
dat <- c(1,0,0,2,1,1,1,2,2,1,1,0,0,1)

## using the constructor
x1 <- new("SNPbin", dat)
x1

## using 'as' methods
x2 <- as.SNPbin(dat)
x3 <- as(dat, "SNPbin")

identical(x1,x2)
identical(x1,x3)
}
}
\keyword{classes}
