\encoding{UTF-8}
\name{makefreq}
\alias{makefreq}
\title{Function to generate allelic frequencies}
\description{The function \code{makefreq} generates a table of allelic
  frequencies from an object of class \code{genpop}.
}
\usage{
makefreq(x,quiet=FALSE,missing=NA,truenames=TRUE)
}
\arguments{
  \item{x}{an object of class \code{genpop}.}
  \item{quiet}{logical stating whether a conversion message must be
    printed (TRUE,default) or not (FALSE).}
  \item{missing}{treatment for missing values. Can be NA, 0 or "mean"
    (see details)}
  \item{truenames}{a logical indicating whether true labels (as opposed
    to generic labels) should be used to name the output.}
}
\details{There are 3 treatments for missing values: \cr
  - NA: kept as NA.\cr
  - 0: missing values are considered as zero. Recommended for a PCA on
  compositionnal data.\cr
  - "mean": missing values are given the mean frequency of the
  corresponding allele. Recommended for a centred PCA.  
}
\value{Returns a list with the following components:
  \item{tab}{matrix of allelic frequencies (rows: populations; columns: alleles).}
  \item{nobs}{number of observations (i.e. alleles) for each population x locus combinaison.}
  \item{call}{the matched call}
}
\seealso{\code{\link{genpop}}
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
data(microbov)
obj1 <- microbov

obj2 <- genind2genpop(obj1)

Xfreq <- makefreq(obj2,missing="mean")


# perform a correspondance analysis on counts data

Xcount <- genind2genpop(obj1,missing="chi2")
ca1 <- dudi.coa(as.data.frame(Xcount@tab),scannf=FALSE)
s.label(ca1$li,sub="Correspondance Analysis",csub=1.2)
add.scatter.eig(ca1$eig,nf=2,xax=1,yax=2,posi="topleft")

# perform a principal component analysis on frequency data
pca1 <- dudi.pca(Xfreq$tab,scale=FALSE,scannf=FALSE)
s.label(pca1$li,sub="Principal Component Analysis",csub=1.2)
add.scatter.eig(pca1$eig,nf=2,xax=1,yax=2,posi="top")
}
}
\keyword{manip}
\keyword{multivariate}
