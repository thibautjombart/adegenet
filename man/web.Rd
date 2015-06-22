\encoding{UTF-8}
\docType{methods}
\name{Website and tutorials}
\alias{adegenetWeb}
\alias{adegenetTutorial}
\alias{adegenetTutorials}
\alias{adegenetIssue}
\alias{adegenetIssues}
\title{ Web resources for adegenet}
\description{
  These functions can be use to access online resources for adegenet.\cr

  They include:
  \itemize{
    \item \code{adegenetWeb}: opens the adegenet website in a web navigator.
    \item \code{adegenetTutorial}: opens a specified tutorial.
    \item \code{adegenetIssues}: opens issues system on github.
  }
}
\usage{
adegenetWeb()

adegenetTutorial(which=c("basics","spca","dapc","genomics"))

adegenetIssues()
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
adegenetIssues()
}
}
