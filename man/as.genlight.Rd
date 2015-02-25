\name{as.genlight}
\alias{as,genlight,matrix-method}
\alias{as,genlight,data.frame-method}
\alias{as,genlight,list-method}
% \alias{as,genlight,snp.matrix-method}
\alias{as.genlight}
\alias{as.genlight,matrix-method}
\alias{as.genlight,data.frame-method}
\alias{as.genlight,list-method}
% \alias{as.genlight,snp.matrix-method}
\alias{coerce,genlight,matrix-method}
\alias{coerce,genlight,data.frame-method}
\alias{coerce,genlight,list-method}
% \alias{coerce,genlight,snp.matrix-method}

% \alias{,genlight-method}
% \alias{,genlight-method}
% \alias{,genlight-method}
% \alias{,genlight-method}
%%%%
\title{Conversion to class "genlight"}
\description{
  The class \code{genlight} is a formal (S4) class for storing a genotypes
  of binary SNPs in a compact way, using a bit-level coding scheme. New
  instances of this class are best created using \code{new}; see the
  manpage of \linkS4class{genlight} for more information on this point.

  As a shortcut, conversion methods can be used to convert various
  objects into a \linkS4class{genlight} object. Conversions can be
  achieved using S3-style (\code{as.genlight(x)}) or S4-style
  (\code{as(x,"genlight"}) procedures. All of them call upon the
  constructor (\code{new}) of \linkS4class{genlight} objects.

  Conversion is currently available from the following objects:
  - matrix of type integer/numeric
  - data.frame with integer/numeric data
  - list of vectors of integer/numeric type
}
\author{Thibaut Jombart (\email{t.jombart@imperial.ac.uk})}
\seealso{
 Related class:\cr
  -  \code{\linkS4class{SNPbin}}, for storing individual genotypes of
  binary SNPs\cr
  
  -  \code{\linkS4class{genind}}
}
\examples{
\dontrun{
## data to be converted
dat <- list(toto=c(1,1,0,0,2,2,1,2,NA), titi=c(NA,1,1,0,1,1,1,0,0), tata=c(NA,0,3, NA,1,1,1,0,0))

## using the constructor
x1 <- new("genlight", dat)
x1

## using 'as' methods
x2 <- as.genlight(dat)
x3 <- as(dat, "genlight")

identical(x1,x2)
identical(x1,x3)
}

}
\keyword{classes}
