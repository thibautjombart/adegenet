\encoding{UTF-8}
\docType{methods}
\name{Website and tutorials}
\alias{adegenetWeb}
\alias{adegenetTutorial}
\alias{adegenetTutorials}
\title{ Website and tutorials for adegenet}
\description{
  These functions allow to access documentation for adegenet available online.\cr

  These functions include:
  \itemize{
    \item \code{adegenetWeb}: opens the adegenet website in a web navigator.
    \item \code{adegenetTutorial}: opens a specified tutorial.
  }
}
\usage{
adegenetWeb()

adegenetTutorial(which=c("basics","spca","dapc","genomics"))
}
\arguments{
  \item{which}{a character string indicating the type of tutorial to open}
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
## this opens the adegenet website
adegenetWeb()
adegenetTutorial("dapc")
}
}
