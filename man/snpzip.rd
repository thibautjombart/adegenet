\encoding{UTF-8}
\name{snpzip}
\alias{snpzip}
\title{Identification of structural SNPs}
\description{
  The function \code{snpzip} identifies the set of alleles which contribute most
  significantly to phenotypic structure.

  This procedure first runs a Discriminant Analysis of Principal Components (DAPC)
  to quantify the contribution of individual alleles to between-population
  structure. Then, defining contribution to DAPC as the measure of distance
  between alleles, hierarchical clustering is used to identify two groups
  of alleles: structural SNPs and non-structural SNPs.
}
\usage{
snpzip(snps,phen,plot=TRUE,pca.plot=FALSE,method=c("complete","single",
"average","centroid","mcquitty","median","ward"), \dots)
}
\arguments{
  \item{snps}{a \code{matrix} used as input of DAPC.}
  \item{phen}{a \code{factor} indicating the group membership of individuals.}
  \item{plot}{a \code{logical} indicating whether a graphical representation of the 
    DAPC results should be displayed.}
  \item{pca.plot}{a \code{logical} indicating whether the results of the 
    cross-validation step should be displayed.}
  \item{method}{the clustering method to be used. This should be 
    (an unambiguous abbreviation of) one of \code{"complete", "single", "average", 
    "centroid", "mcquitty", "median",} or \code{"ward"}.} 
  \item{\dots}{further arguments.}
    %not sure I have any further arguments?
}

\details{
  \code{snpzip} provides an objective procedure to delineate between structural 
  and non-structural SNPs identified by Discriminant Analysis of Principal Components 
  (DAPC, Jombart et al. 2010). 
  \code{snpzip} precedes the multivariate analysis with a cross-validation step 
  to ensure that the subsequent DAPC is performed optimally.
  The contributions of alleles to the DAPC are then submitted to \code{hclust}, 
  where they define a distance matrix upon which hierarchical clustering is carried out.
  To complete the procedure, \code{snpzip} uses \code{cutree} to automatically 
  subdivide the set of SNPs fed into the analysis into two groups: 
  those which contribute significantly to the phenotypic structure of interest, 
  and those which do not.   
}

\value{
  A \code{list} with four items: 
  the first cites the number of principal components (PCs) of PCA retained in the DAPC, 
  the second indicates the number of structural and non-structural SNPs identified by 
  \code{snpzip}, the third provides a list of the structuring alleles, and the 
  fourth item details the contributions of these structuring alleles to the DAPC.
  
  If \code{plot=TRUE}, a scatter plot will provide a visualization of the DAPC results.
  
  If \code{pca.plot=TRUE}, the results of the cross-validation step will be displayed 
  as an \code{array} showing mean assignment success by number of PCs of PCA; 
  additionally, the number of PCs producing the highest mean success will be indicated, 
  the distribution of mean sucess for random choice will be presented in 
  an \code{array} of quantiles, and a scatter plot of the results of 
  cross-validation will be provided.   
}

\references{
Jombart T, Devillard S and Balloux F (2010) Discriminant analysis of principal 
components: a new method for the analysis of genetically structured populations. 
BMC Genetics11:94. doi:10.1186/1471-2156-11-94
}

\author{ Caitlin Collins \email{caitlin.collins12@imperial.ac.uk} }
\examples{
\dontrun{
## Not sure what examples I can usefully provide here. 
}
}
\keyword{multivariate}