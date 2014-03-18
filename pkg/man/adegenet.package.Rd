\encoding{UTF-8}
\name{adegenet-package}
\alias{adegenet-package}
\alias{adegenet}
\docType{package}
\title{The adegenet package}

\description{

  This package is devoted to the multivariate analysis of
  genetic markers data. These data can be codominant markers (e.g. microsatellites) or
  presence/absence data (e.g. AFLP), and have any level of ploidy.
  'adegenet' defines three formal (S4) classes:\cr
  - \linkS4class{genind}: a class for data of individuals ("genind" stands for genotypes-individuals).\cr
  - \linkS4class{genpop}: a class for data of groups of individuals
  ("genpop" stands for genotypes-populations)\cr
  - \linkS4class{genlight}: a class for genome-wide SNP data\cr

  For more information about these classes, type "class ? genind",
  "class ? genpop", or "?genlight".\cr

  Essential functionalities of the package are presented througout 4
  tutorial vignettes, accessible using \code{vignette("name-below", package="adegenet")}:\cr
  - adegenet-basics: introduction to the package.\cr
  - adegenet-spca: multivariate analysis of spatial genetic patterns.\cr
  - adegenet-dapc: population structure and group assignment using DAPC.\cr
  - adegenet-genomics: introduction to the class \linkS4class{genlight}
  for the handling and analysis of genome-wide SNP data.\cr

  Important functions are also summarized below.\cr

  
  === IMPORTING DATA ===\cr
  = TO GENIND OBJECTS = \cr
  \code{adegenet} imports data to \linkS4class{genind} object from the
  following softwares:\cr
  - STRUCTURE: see \code{\link{read.structure}}\cr
  - GENETIX: see \code{\link{read.genetix}}\cr
  - FSTAT: see \code{\link{read.fstat}}\cr
  - Genepop: see \code{\link{read.genepop}}\cr
  To import data from any of these formats, you can also use the general
  function \code{\link{import2genind}}.\cr

  In addition, it can extract polymorphic sites from nucleotide and
  amino-acid alignments:\cr
  - DNA files: use  \code{\link[ape]{read.dna}} from the ape package,
  and then extract SNPs from DNA alignments using
  \code{\link{DNAbin2genind}}. \cr

 - protein sequences alignments: polymorphic sites can be extracted from
   protein sequences alignments in \code{alignment} format (package
   \code{seqinr}, see \code{\link[seqinr]{as.alignment}}) using the
   function \code{\link{alignment2genind}}. \cr

   The function \code{\link{fasta2DNAbin}} allows for reading fasta
   files into DNAbin object with minimum RAM requirements.\cr

  It is also possible to read genotypes coded by character strings from
  a data.frame in which genotypes are in rows, markers in columns. For
  this, use \code{\link{df2genind}}. Note that \code{\link{df2genind}}
  can be used for any level of ploidy.\cr

  = TO GENLIGHT OBJECTS = \cr
  SNP data can be read from the following formats:\cr
  - PLINK: see function \code{\link{read.PLINK}}\cr
  - .snp (adegenet's own format): see function \code{\link{read.snp}}\cr

  SNP can also be extracted from aligned DNA sequences with the fasta
  format, using \code{\link{fasta2genlight}}\cr


   === EXPORTING DATA ===\cr
   \code{adegenet} exports data from \linkS4class{genind} object to
   formats recognized by other R packages:\cr
   - the genetics package: see \code{\link{genind2genotype}}\cr
   - the hierfstat package: see \code{\link{genind2hierfstat}}\cr

   Genotypes can also be recoded from a \linkS4class{genind} object into
   a data.frame of character strings, using any separator between
   alleles. This covers formats from many softwares like GENETIX or
   STRUCTURE. For this, see \code{\link{genind2df}}.\cr

   Also note that the \code{pegas} package imports \linkS4class{genind}
   objects using the function \code{as.loci}.

   
   === MANIPULATING DATA ===\cr
   Several functions allow one to manipulate \linkS4class{genind} or
   \linkS4class{genpop} objects\cr
   - \code{\link{genind2genpop}}: convert a \linkS4class{genind} object
   to a  \linkS4class{genpop} \cr
   - \code{\link{seploc}}: creates one object per marker; for
   \linkS4class{genlight} objects, creates blocks of SNPs.\cr
   - \code{\link{seppop}}: creates one object per population \cr
   - \code{\link{na.replace}}: replaces missing data (NA) in an
   approriate way \cr
   - \code{\link{truenames}}: restores true names of an object
   (\linkS4class{genind} and \linkS4class{genpop} use generic labels)
   \cr
   - x[i,j]: create a new object keeping only genotypes (or populations)
   indexed by 'i' and the alleles indexed by 'j'.\cr
   - \code{\link{makefreq}}: returns a table of allelic frequencies from
   a \linkS4class{genpop} object.\cr
   - \code{\link{repool}} merges genoptypes from different
   gene pools into one single \linkS4class{genind} object.\cr
   - \code{\link{propTyped}} returns the proportion of available (typed)
   data, by individual, population, and/or locus.\cr
   - \code{\link{selPopSize}} subsets data, retaining only genotypes
   from a population whose sample size is above a given level.\cr
   - \code{\link{pop}} sets the population of a set of genotypes.\cr


   === ANALYZING DATA ===\cr
   Several functions allow to use usual, and less usual analyses:\cr
   - \code{\link{HWE.test.genind}}: performs HWE test for all
   populations and loci combinations \cr
   - \code{\link{pairwise.fst}}: computes simple pairwise Fst between populations\cr
   - \code{\link{dist.genpop}}: computes 5 genetic distances among populations. \cr
   - \code{\link{monmonier}}: implementation of the Monmonier algorithm,
   used to seek genetic boundaries among individuals or
   populations. Optimized boundaries can be obtained using
   \code{\link{optimize.monmonier}}. Object of the class
   \code{monmonier} can be plotted and printed using the corresponding
   methods. \cr
   - \code{\link{spca}}: implements Jombart et al. (2008) spatial
   Principal Component Analysis \cr
   - \code{\link{global.rtest}}: implements Jombart et al. (2008)
   test for global spatial structures \cr
   - \code{\link{local.rtest}}: implements Jombart et al. (2008)
   test for local spatial structures \cr
   - \code{\link{propShared}}: computes the proportion of shared
   alleles in a set of genotypes (i.e. from a genind object)\cr
   - \code{\link{propTyped}}: function to investigate missing data in
   several ways \cr
   - \code{\link{scaleGen}}: generic method to scale
   \linkS4class{genind} or \linkS4class{genpop} before a principal
   component analysis \cr
   - \code{\link{Hs}}: computes the average expected heterozygosity by
   population in a \linkS4class{genpop}. Classically Used as a measure
   of genetic diversity.\cr
   - \code{\link{find.clusters}} and \code{\link{dapc}}: implement the
   Discriminant Analysis of Principal Component (DAPC, Jombart et al.,
   2010).\cr
   - \code{\link{seqTrack}}: implements the SeqTrack algorithm for
   recontructing transmission trees of pathogens (Jombart et al.,
   2010) .\cr
   \code{\link{glPca}}: implements PCA for \linkS4class{genlight}
   objects.\cr
   - \code{\link{gengraph}}: implements some simple graph-based
   clustering using genetic data.
   - \code{\link{snpposi.plot}} and \code{\link{snpposi.test}}:
   visualize the distribution of SNPs on a genetic sequence and test
   their randomness.
   - \code{\link{adegenetServer}}: opens up a web interface for some
   functionalities of the package (DAPC with cross validation and
   feature selection).\cr



   === GRAPHICS ===\cr
   - \code{\link{colorplot}}: plots points with associated values for up
   to three variables represented by colors using the RGB system;
   useful for spatial mapping of principal components.\cr
   - \code{\link{loadingplot}}: plots loadings of variables. Useful for
   representing the contribution of alleles to a given principal
   component in a multivariate method. \cr
 - \code{\link{scatter.dapc}}: scatterplots for DAPC results.\cr 
   - \code{\link{compoplot}}: plots membership probabilities from a DAPC
   object. \cr


   === SIMULATING DATA ===\cr
   - \code{\link{hybridize}}: implements hybridization between two populations. \cr
   - \code{\link{haploGen}}: simulates genealogies of haplotypes,
   storing full genomes. \cr
   % - \code{\link{haploPop}}: simulates populations of haplotypes, using
   % different population dynamics, storing SNPs (under development). \cr
   - \code{\link{glSim}}: simulates simple \linkS4class{genlight} objects.\cr


   === DATASETS ===\cr
   - \code{\link{H3N2}}: Seasonal influenza (H3N2) HA segment data. \cr
   - \code{\link{dapcIllus}}: Simulated data illustrating the DAPC. \cr
   - \code{\link{eHGDP}}: Extended HGDP-CEPH dataset. \cr
   - \code{\link{microbov}}: Microsatellites genotypes of 15 cattle breeds. \cr
   - \code{\link{nancycats}}: Microsatellites genotypes of 237 cats from 17 colonies of Nancy (France). \cr
   - \code{\link{rupica}}: Microsatellites genotypes of 335 chamois
   (Rupicapra rupicapra) from the Bauges mountains (France).\cr
   - \code{\link{sim2pop}}: Simulated genotypes of two georeferenced populations.\cr
   - \code{\link{spcaIllus}}: Simulated data illustrating the sPCA. \cr

   For more information, visit the adegenet website using the function
   \code{\link{adegenetWeb}}.\cr

   Tutorials are available via the command \code{\link{adegenetTutorials}}.\cr


   To cite adegenet, please use the reference given by
   \code{citation("adegenet")} (or see reference below).
 }
\details{
  \tabular{ll}{
    Package: \tab adegenet\cr
    Type: \tab Package\cr
    Version: \tab 1.4-1\cr
    Date: \tab 2014-03-14 \cr
    License: \tab GPL (>=2)
  }
}
\author{
  Thibaut Jombart <t.jombart@imperial.ac.uk>\cr
  Developers: Caitlin Collins <caitiecollins17@gmail.com> Ismail Ahmed
  <ismail.ahmed@inserm.fr>, Federico Calboli <f.calboli@imperial.ac.uk>,
  Tobias Erik Reiners, Peter Solymos, Anne Cori\cr
  Contributed datasets from: Katayoun Moazami-Goudarzi, Denis Laloë,
  Dominique Pontier, Daniel Maillard, Francois Balloux.  }
\references{
  Jombart T. (2008) adegenet: a R package for the multivariate analysis
  of genetic markers \emph{Bioinformatics} 24: 1403-1405. doi:
  10.1093/bioinformatics/btn129\cr

  Jombart T. and Ahmed I. (2011) adegenet 1.3-1: new tools for the
  analysis of genome-wide SNP data.  \emph{Bioinformatics}. doi:
  10.1093/bioinformatics/btr521

  Jombart T, Devillard S and Balloux F (2010) Discriminant analysis of
  principal components: a new method for the analysis of genetically
  structured populations. BMC Genetics 11:94.
  doi:10.1186/1471-2156-11-94\cr

  Jombart T, Eggo R, Dodd P, Balloux F (2010) Reconstructing disease
  outbreaks from genetic data: a graph approach. \emph{Heredity}. doi:
  10.1038/hdy.2010.78.\cr

  Jombart, T., Devillard, S., Dufour, A.-B. and Pontier, D. (2008) Revealing
  cryptic spatial patterns in genetic variability by a new multivariate
  method. \emph{Heredity}, \bold{101}, 92--103.\cr

  See adegenet website: \url{http://adegenet.r-forge.r-project.org/}\cr

  Please post your questions on 'the adegenet forum': adegenet-forum@lists.r-forge.r-project.org

}
\keyword{manip}
\keyword{multivariate}
\seealso{
  adegenet is related to several packages, in particular:\cr
  - \code{ade4} for multivariate analysis\cr
  - \code{pegas} for population genetics tools\cr
  - \code{ape} for phylogenetics and DNA data handling\cr
  - \code{seqinr} for handling nucleic and proteic sequences\cr
  - \code{shiny} for R-based web interfaces\cr
}

