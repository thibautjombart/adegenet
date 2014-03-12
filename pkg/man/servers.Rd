\encoding{UTF-8}
\name{Adegenet servers}
\alias{adegenetServer}
\title{ Web servers for adegenet }
\description{
  The function \code{adegenetServer} opens up a web page providing a
  simple user interface for some of the functionalities implemented in
  adegenet. These servers have been developed using the package
  \code{shiny}.\cr

  Currently available servers include:
  \itemize{
    \item \code{DAPC}: a server for the Discriminant Analysis of
    Principal Components (see ?dapc)
  }
}
\seealso{
  \link{dapc}
}
\usage{
adegenetServer(what=c("DAPC"))
}
\arguments{
  \item{what}{a character string indicating which server to start;
  currently accepted values are: "DAPC"}
}
\value{
 The function invisibly returns NULL.
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}
  Caitlin Collins
}
\examples{
\dontrun{
## this opens a web page for DAPC
adegenetServer()
}
}
