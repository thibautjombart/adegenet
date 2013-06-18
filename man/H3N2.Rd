\encoding{UTF-8}
\name{H3N2}
\alias{H3N2}
\alias{usflu}
\alias{usflu.fasta}
\alias{USflu}
\alias{USflu.fasta}
\docType{data}
\title{Seasonal influenza (H3N2) HA segment data}
\description{
The dataset \code{H3N2} consists of 1903 strains of seasonal influenza
(H3N2) distributed worldwide, and typed at 125 SNPs located in the
hemagglutinin (HA) segment. It is stored as an R object with class
\linkS4class{genind} and can be accessed as usual using
\code{data(H3N2)} (see example). These data were gathered from DNA
sequences available from Genbank (http://www.ncbi.nlm.nih.gov/Genbank/).


The data file \code{usflu.fasta} is a toy dataset also gathered from
Genbank, consisting of the aligned sequences of 80 seasonal influenza
isolates (HA segment) sampled in the US, in \code{fasta} format. This
file is installed alongside the package; the path to this file is
automatically determined by R using \code{system.file} (see example in
this manpage and in ?fasta2genlight) as well.
}
\usage{data(H3N2)}
\format{
    \code{H3N2} is a genind object with several data frame as supplementary
    components (\code{H3N2@other) slort}, which contains the following items:
    \describe{
        \item{x}{a \code{data.frame} containing miscellanous annotations
	  of the sequences.}
        \item{xy}{a matrix with two columns indicating the geographic
	  coordinates of the strains, as longitudes and latitudes.}
        \item{epid}{a character vector indicating the epidemic of the strains.}
    }
}
\source{
  This dataset was prepared by Thibaut Jombart
  (t.jombart@imperia.ac.uk), from annotated sequences available on
  Genbank (http://www.ncbi.nlm.nih.gov/Genbank/).
}
\references{
Jombart, T., Devillard, S. and Balloux, F.
Discriminant analysis of principal components: a new method for the analysis of
genetically structured populations. Submitted to \emph{BMC genetics}.
}
\examples{
\dontrun{
#### H3N2 ####
## LOAD DATA
data(H3N2)
H3N2

## set population to yearly epidemics
pop(H3N2) <- factor(H3N2$other$epid)



## PERFORM DAPC - USE POPULATIONS AS CLUSTERS
## to reproduce exactly analyses from the paper, use "n.pca=1000"
dapc1 <- dapc(H3N2, all.contrib=TRUE, scale=FALSE, n.pca=150, n.da=5)
dapc1

## (see ?dapc for details about the output)


## SCREEPLOT OF EIGENVALUES
barplot(dapc1$eig, main="H3N2 - DAPC eigenvalues")


## SCATTERPLOT (axes 1-2)
scatter(dapc1, posi.da="topleft", cstar=FALSE, cex=2, pch=17:22,
solid=.5, bg="white")




#### usflu.fasta ####
myPath <- system.file("files/usflu.fasta",package="adegenet")
myPath

## extract SNPs from alignments using fasta2genlight
## see ?fasta2genlight for more details
obj <- fasta2genlight(myPath, chunk=10) # process 10 sequences at a time
obj
}
}
\keyword{datasets}
