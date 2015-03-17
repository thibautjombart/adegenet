

#' Accessors for adegenet objects
#' 
#' An accessor is a function that allows to interact with slots of an object in
#' a convenient way. Several accessors are available for \linkS4class{genind}
#' or \linkS4class{genpop} objects. The operator "\$" and "\$<-" are used to
#' access the slots, being equivalent to "@@" and "@@<-".\cr
#' 
#' The operator "[" can be used to access components of the matrix slot "@@tab",
#' returning a \linkS4class{genind} or \linkS4class{genpop} object. This syntax
#' is the same as for a matrix; for instance:\cr - "obj[,]" returns "obj" \cr -
#' "obj[1:10,]" returns an object with only the first 10 genotypes (if "obj" is
#' a \linkS4class{genind}) or the first 10 populations (if "obj" is a
#' \linkS4class{genpop}) of "obj" \cr - "obj[1:10, 5:10]" returns an object
#' keeping the first 10 entities and the alleles 5 to 10.\cr -
#' "obj[loc=c("L1","L3")]" returns an object keeping only the loci specified in
#' the \code{loc} argument (using generic names, not true names; in this
#' example, only the first and the third locus would be retained)\cr -
#' "obj[1:3, drop=TRUE]" returns the first 3 genotypes/populations of "obj",
#' but retaining only alleles that are present in this subset (as opposed to
#' keeping all alleles of "obj", which is the default behavior).\cr
#' 
#' The argument \code{treatOther} handles the treatment of objects in the
#' \code{@@other} slot (see details). The argument \code{drop} can be set to
#' TRUE to drop alleles that are no longer represented in the subset.
#' 
#' The "[" operator can treat elements in the \code{@@other} slot as well. For
#' instance, if \code{obj@@other$xy} contains spatial coordinates, the
#' \code{obj[1:3,]@@other$xy} will contain the spatial coordinates of the
#' genotypes (or population) 1,2 and 3. This is handled through the argument
#' \code{treatOther}, a logical defaulting to TRUE. If set to FALSE, the
#' \code{@@other} returned unmodified.\cr
#' 
#' Note that only matrix-like, vector-like and lists can be proceeded in
#' \code{@@other}. Other kind of objects will issue a warning an be returned as
#' they are, unless the argument \code{quiet} is left to TRUE, its default
#' value.\cr
#' 
#' The \code{drop} argument can be set to TRUE to retain only alleles that are
#' present in the subset. To achieve better control of polymorphism of the
#' data, see \code{\link{isPoly}}.
#' 
#' @name Accessors
#' @aliases $,genind-method $,genpop-method $<-,genind-method $<-,genpop-method
#' [,genind-method [,genpop-method nLoc nLoc,genind-method nLoc,genpop-method
#' nInd nInd,genind-method pop pop<- pop,genind-method pop<-,genind-method
#' locNames locNames,genind-method locNames,genpop-method locNames<-
#' locNames<-,genind-method locNames<-,genpop-method indNames
#' indNames,genind-method indNames<- indNames<-,genind-method ploidy
#' ploidy,genind-method ploidy,genpop-method ploidy<- ploidy<-,genind-method
#' ploidy<-,genpop-method alleles alleles,genind-method alleles,genpop-method
#' alleles<- alleles<-,genind-method alleles<-,genpop-method other
#' other,genind-method other,genpop-method other<- other<-,genind-method
#' other<-,genpop-method
#' @docType methods
#' @param x a \linkS4class{genind} or a \linkS4class{genpop} object.
#' @param withAlleles a logical indicating whether the result should be of the
#' form [locus name].[allele name], instead of [locus name].
#' @param \dots further arguments to be passed to other methods (currently not
#' used).
#' @return A \linkS4class{genind} or \linkS4class{genpop} object.
#' @section Methods: \describe{ \item{nInd}{returns the number of individuals
#' in the \code{genind} object} \item{nLoc}{returns the number of loci of the
#' object} \item{pop}{returns the population factor of the object, using true
#' (as opposed to generic) levels.} \item{pop<-}{replacement method for the
#' \code{@@pop} slot of an object. The content of \code{@@pop} and
#' \code{@@pop.names} is updated automatically.} \item{indNames}{returns the
#' true names of individuals.} \item{indNames<-}{sets the true names of
#' individuals using a vector of length \code{nInd(x)}.} \item{locNames}{returns the
#' true names of markers and/or alleles.} \item{locNames<-}{sets the true names
#' of markers using a vector of length \code{nLoc(x)}.} \item{ploidy}{returns the
#' ploidy of the data.} \item{ploidy<-}{sets the ploidy of the data using an
#' integer.} \item{alleles}{returns the alleles of each locus.}
#' \item{alleles<-}{sets the alleles of each locus using a list with one
#' character vector for each locus.} \item{other}{returns the content of the
#' \code{@@other} slot (misc. information); returns \code{NULL} if the slot is
#' empty or of length zero.} \item{other<-}{sets the content of the
#' \code{@@other} slot (misc. information); the provided value needs to be a
#' list; it not, provided value will be stored within a list.} }
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords manip
#' @examples
#' 
#' data(nancycats)
#' nancycats
#' pop(nancycats) # get the populations
#' indNames(nancycats) # get the labels of individuals
#' locNames(nancycats) # get the labels of the loci
#' alleles(nancycats) # get the alleles
#' 
#' # let's isolate populations 4 and 8
#' temp <- nancycats@@pop=="P04" | nancycats@@pop=="P08"
#' obj <- nancycats[temp,]
#' obj
#' 
#' pop(obj)
#' 
#' # let's isolate two markers, fca23 and fca90
#' locNames(nancycats)
#' 
#' # they correspond to L2 and L7
#' nancycats$loc.fac
#' temp <- nancycats$loc.fac=="L2" | nancycats$loc.fac=="L7"
#' obj <- nancycats[,temp]
#' obj
#' 
#' obj$loc.fac
#' locNames(obj)
#' 
#' # or more simply
#' nancycats[loc=c("L2","L7")]
#' obj$loc.fac 
#' locNames(obj)
#' 
#' # using 'drop':
#' truenames(nancycats[1:2])$tab
#' truenames(nancycats[1:2, drop=TRUE])$tab
#' 
#' # illustrate how 'other' slot is handled
#' colonies <- genind2genpop(nancycats)
#' colonies@@other$aChar <- "This will not be proceeded"
#' colonies123 <- colonies[1:3]
#' colonies
#' colonies@@other$xy
#' 
#' # illustrate pop
#' obj <- nancycats[sample(1:100,10)]
#' obj$pop
#' obj$pop.names
#' pop(obj)
#' pop(obj) <- rep(c('b','a'), each=5)
#' obj$pop
#' obj$pop.names
#' pop(obj)
#' 
#' # illustrate locNames
#' locNames(obj)
#' locNames(obj, withAlleles=TRUE)
#' 
#' 
NULL





#' The adegenet package
#' 
#' 
#' This package is devoted to the multivariate analysis of genetic markers
#' data. These data can be codominant markers (e.g. microsatellites) or
#' presence/absence data (e.g. AFLP), and have any level of ploidy.  'adegenet'
#' defines three formal (S4) classes:\cr - \linkS4class{genind}: a class for
#' data of individuals ("genind" stands for genotypes-individuals).\cr -
#' \linkS4class{genpop}: a class for data of groups of individuals ("genpop"
#' stands for genotypes-populations)\cr - \linkS4class{genlight}: a class for
#' genome-wide SNP data\cr
#' 
#' For more information about these classes, type "class ? genind", "class ?
#' genpop", or "?genlight".\cr
#' 
#' Essential functionalities of the package are presented througout 4
#' tutorials, accessible using \code{adegenetTutorial(which="name-below")}:\cr
#' - \code{basics}: introduction to the package.\cr - \code{spca}: multivariate
#' analysis of spatial genetic patterns.\cr - \code{dapc}: population structure
#' and group assignment using DAPC.\cr - \code{genomics}: introduction to the
#' class \linkS4class{genlight} for the handling and analysis of genome-wide
#' SNP data.\cr
#' 
#' Note: In older versions of adegenet, these tutorials were avilable as
#' vignettes, accessible through the function \code{vignette("name-below",
#' package="adegenet")}:\cr - \code{adegenet-basics}.\cr -
#' \code{adegenet-spca}.\cr - \code{adegenet-dapc}.\cr -
#' \code{adegenet-genomics}.\cr
#' 
#' Important functions are also summarized below.\cr
#' 
#' === IMPORTING DATA ===\cr = TO GENIND OBJECTS = \cr \code{adegenet} imports
#' data to \linkS4class{genind} object from the following softwares:\cr -
#' STRUCTURE: see \code{\link{read.structure}}\cr - GENETIX: see
#' \code{\link{read.genetix}}\cr - FSTAT: see \code{\link{read.fstat}}\cr -
#' Genepop: see \code{\link{read.genepop}}\cr To import data from any of these
#' formats, you can also use the general function
#' \code{\link{import2genind}}.\cr
#' 
#' In addition, it can extract polymorphic sites from nucleotide and amino-acid
#' alignments:\cr - DNA files: use \code{\link[ape]{read.dna}} from the ape
#' package, and then extract SNPs from DNA alignments using
#' \code{\link{DNAbin2genind}}. \cr
#' 
#' - protein sequences alignments: polymorphic sites can be extracted from
#' protein sequences alignments in \code{alignment} format (package
#' \code{seqinr}, see \code{\link[seqinr]{as.alignment}}) using the function
#' \code{\link{alignment2genind}}. \cr
#' 
#' The function \code{\link{fasta2DNAbin}} allows for reading fasta files into
#' DNAbin object with minimum RAM requirements.\cr
#' 
#' It is also possible to read genotypes coded by character strings from a
#' data.frame in which genotypes are in rows, markers in columns. For this, use
#' \code{\link{df2genind}}. Note that \code{\link{df2genind}} can be used for
#' any level of ploidy.\cr
#' 
#' = TO GENLIGHT OBJECTS = \cr SNP data can be read from the following
#' formats:\cr - PLINK: see function \code{\link{read.PLINK}}\cr - .snp
#' (adegenet's own format): see function \code{\link{read.snp}}\cr
#' 
#' SNP can also be extracted from aligned DNA sequences with the fasta format,
#' using \code{\link{fasta2genlight}}\cr
#' 
#' === EXPORTING DATA ===\cr \code{adegenet} exports data from
#' \linkS4class{genind} object to formats recognized by other R packages:\cr -
#' the genetics package: see \code{\link{genind2genotype}}\cr - the hierfstat
#' package: see \code{\link{genind2hierfstat}}\cr
#' 
#' Genotypes can also be recoded from a \linkS4class{genind} object into a
#' data.frame of character strings, using any separator between alleles. This
#' covers formats from many softwares like GENETIX or STRUCTURE. For this, see
#' \code{\link{genind2df}}.\cr
#' 
#' Also note that the \code{pegas} package imports \linkS4class{genind} objects
#' using the function \code{as.loci}.
#' 
#' === MANIPULATING DATA ===\cr Several functions allow one to manipulate
#' \linkS4class{genind} or \linkS4class{genpop} objects\cr -
#' \code{\link{genind2genpop}}: convert a \linkS4class{genind} object to a
#' \linkS4class{genpop} \cr - \code{\link{seploc}}: creates one object per
#' marker; for \linkS4class{genlight} objects, creates blocks of SNPs.\cr -
#' \code{\link{seppop}}: creates one object per population \cr -
#' \code{\link{na.replace}}: replaces missing data (NA) in an approriate way
#' \cr - \code{\link{truenames}}: restores true names of an object
#' (\linkS4class{genind} and \linkS4class{genpop} use generic labels) \cr -
#' x[i,j]: create a new object keeping only genotypes (or populations) indexed
#' by 'i' and the alleles indexed by 'j'.\cr - \code{\link{makefreq}}: returns
#' a table of allelic frequencies from a \linkS4class{genpop} object.\cr -
#' \code{\link{repool}} merges genoptypes from different gene pools into one
#' single \linkS4class{genind} object.\cr - \code{\link{propTyped}} returns the
#' proportion of available (typed) data, by individual, population, and/or
#' locus.\cr - \code{\link{selPopSize}} subsets data, retaining only genotypes
#' from a population whose sample size is above a given level.\cr -
#' \code{\link{pop}} sets the population of a set of genotypes.\cr
#' 
#' === ANALYZING DATA ===\cr Several functions allow to use usual, and less
#' usual analyses:\cr - \code{\link{HWE.test.genind}}: performs HWE test for
#' all populations and loci combinations \cr - \code{\link{pairwise.fst}}:
#' computes simple pairwise Fst between populations\cr -
#' \code{\link{dist.genpop}}: computes 5 genetic distances among populations.
#' \cr - \code{\link{monmonier}}: implementation of the Monmonier algorithm,
#' used to seek genetic boundaries among individuals or populations. Optimized
#' boundaries can be obtained using \code{\link{optimize.monmonier}}. Object of
#' the class \code{monmonier} can be plotted and printed using the
#' corresponding methods. \cr - \code{\link{spca}}: implements Jombart et al.
#' (2008) spatial Principal Component Analysis \cr -
#' \code{\link{global.rtest}}: implements Jombart et al. (2008) test for global
#' spatial structures \cr - \code{\link{local.rtest}}: implements Jombart et
#' al. (2008) test for local spatial structures \cr - \code{\link{propShared}}:
#' computes the proportion of shared alleles in a set of genotypes (i.e. from a
#' genind object)\cr - \code{\link{propTyped}}: function to investigate missing
#' data in several ways \cr - \code{\link{scaleGen}}: generic method to scale
#' \linkS4class{genind} or \linkS4class{genpop} before a principal component
#' analysis \cr - \code{\link{Hs}}: computes the average expected
#' heterozygosity by population in a \linkS4class{genpop}. Classically Used as
#' a measure of genetic diversity.\cr - \code{\link{find.clusters}} and
#' \code{\link{dapc}}: implement the Discriminant Analysis of Principal
#' Component (DAPC, Jombart et al., 2010).\cr - \code{\link{seqTrack}}:
#' implements the SeqTrack algorithm for recontructing transmission trees of
#' pathogens (Jombart et al., 2010) .\cr \code{\link{glPca}}: implements PCA
#' for \linkS4class{genlight} objects.\cr - \code{\link{gengraph}}: implements
#' some simple graph-based clustering using genetic data.  -
#' \code{\link{snpposi.plot}} and \code{\link{snpposi.test}}: visualize the
#' distribution of SNPs on a genetic sequence and test their randomness.  -
#' \code{\link{adegenetServer}}: opens up a web interface for some
#' functionalities of the package (DAPC with cross validation and feature
#' selection).\cr
#' 
#' === GRAPHICS ===\cr - \code{\link{colorplot}}: plots points with associated
#' values for up to three variables represented by colors using the RGB system;
#' useful for spatial mapping of principal components.\cr -
#' \code{\link{loadingplot}}: plots loadings of variables. Useful for
#' representing the contribution of alleles to a given principal component in a
#' multivariate method. \cr - \code{\link{scatter.dapc}}: scatterplots for DAPC
#' results.\cr - \code{\link{compoplot}}: plots membership probabilities from a
#' DAPC object. \cr
#' 
#' === SIMULATING DATA ===\cr - \code{\link{hybridize}}: implements
#' hybridization between two populations. \cr - \code{\link{haploGen}}:
#' simulates genealogies of haplotypes, storing full genomes. \cr % -
#' \code{\link{haploPop}}: simulates populations of haplotypes, using %
#' different population dynamics, storing SNPs (under development). \cr -
#' \code{\link{glSim}}: simulates simple \linkS4class{genlight} objects.\cr
#' 
#' === DATASETS ===\cr - \code{\link{H3N2}}: Seasonal influenza (H3N2) HA
#' segment data. \cr - \code{\link{dapcIllus}}: Simulated data illustrating the
#' DAPC. \cr - \code{\link{eHGDP}}: Extended HGDP-CEPH dataset. \cr -
#' \code{\link{microbov}}: Microsatellites genotypes of 15 cattle breeds. \cr -
#' \code{\link{nancycats}}: Microsatellites genotypes of 237 cats from 17
#' colonies of Nancy (France). \cr - \code{\link{rupica}}: Microsatellites
#' genotypes of 335 chamois (Rupicapra rupicapra) from the Bauges mountains
#' (France).\cr - \code{\link{sim2pop}}: Simulated genotypes of two
#' georeferenced populations.\cr - \code{\link{spcaIllus}}: Simulated data
#' illustrating the sPCA. \cr
#' 
#' For more information, visit the adegenet website using the function
#' \code{\link{adegenetWeb}}.\cr
#' 
#' Tutorials are available via the command \code{\link{adegenetTutorials}}.\cr
#' 
#' To cite adegenet, please use the reference given by
#' \code{citation("adegenet")} (or see reference below).
#' 
#' \tabular{ll}{ Package: \tab adegenet\cr Type: \tab Package\cr Version: \tab
#' 1.4-2\cr Date: \tab 2014-05-13 \cr License: \tab GPL (>=2) }
#' 
#' @name adegenet-package
#' @aliases adegenet-package adegenet
#' @docType package
#' @author Thibaut Jombart <t.jombart@@imperial.ac.uk>\cr Developers: Caitlin
#' Collins <caitiecollins17@@gmail.com> Ismail Ahmed <ismail.ahmed@@inserm.fr>,
#' Federico Calboli <f.calboli@@imperial.ac.uk>, Tobias Erik Reiners, Peter
#' Solymos, Anne Cori, Zhian N. Kamvar\cr Contributed datasets from: Katayoun
#' Moazami-Goudarzi, Denis Laloë, Dominique Pontier, Daniel Maillard, Francois
#' Balloux.
#' @seealso adegenet is related to several packages, in particular:\cr -
#' \code{ade4} for multivariate analysis\cr - \code{pegas} for population
#' genetics tools\cr - \code{ape} for phylogenetics and DNA data handling\cr -
#' \code{seqinr} for handling nucleic and proteic sequences\cr - \code{shiny}
#' for R-based web interfaces\cr
#' @references Jombart T. (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers \emph{Bioinformatics} 24: 1403-1405. doi:
#' 10.1093/bioinformatics/btn129\cr
#' 
#' Jombart T. and Ahmed I. (2011) adegenet 1.3-1: new tools for the analysis of
#' genome-wide SNP data.  \emph{Bioinformatics}. doi:
#' 10.1093/bioinformatics/btr521
#' 
#' Jombart T, Devillard S and Balloux F (2010) Discriminant analysis of
#' principal components: a new method for the analysis of genetically
#' structured populations. BMC Genetics 11:94.  doi:10.1186/1471-2156-11-94\cr
#' 
#' Jombart T, Eggo R, Dodd P, Balloux F (2010) Reconstructing disease outbreaks
#' from genetic data: a graph approach. \emph{Heredity}. doi:
#' 10.1038/hdy.2010.78.\cr
#' 
#' Jombart, T., Devillard, S., Dufour, A.-B. and Pontier, D. (2008) Revealing
#' cryptic spatial patterns in genetic variability by a new multivariate
#' method. \emph{Heredity}, \bold{101}, 92--103.\cr
#' 
#' See adegenet website: \url{http://adegenet.r-forge.r-project.org/}\cr
#' 
#' Please post your questions on 'the adegenet forum':
#' adegenet-forum@@lists.r-forge.r-project.org
#' @keywords manip multivariate
NULL





#' Converting genind/genpop objects to other classes
#' 
#' These S3 and S4 methods are used to coerce \linkS4class{genind} and
#' \linkS4class{genpop} objects to matrix-like objects. In most cases, this is
#' equivalent to calling the \code{@@tab} slot. An exception to this is the
#' convertion to \code{\link[ade4]{ktab}} objects used in the ade4 package as
#' inputs for K-tables methods (e.g. Multiple Coinertia Analysis).\cr
#' 
#' 
#' @name as methods in adegenet
#' @aliases as-method as,genind,data.frame-method as,genpop,data.frame-method
#' as,genind,matrix-method as,genpop,matrix-method as,genind,genpop-method
#' ktab-class as,genind,ktab-method as,genpop,ktab-method
#' coerce,genind,data.frame-method coerce,genpop,data.frame-method
#' coerce,genind,matrix-method coerce,genpop,matrix-method
#' coerce,genind,genpop-method coerce,genind,ktab-method
#' coerce,genpop,ktab-method as.data.frame.genind as.data.frame.genpop
#' as.matrix.genind as.matrix.genpop as.genpop.genind as.ktab.genind
#' as.ktab.genpop
#' @docType methods
#' @section Usage: \code{as(object, Class)}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords methods
#' @examples
#' 
#' \dontrun{
#' data(microbov)
#' x <- na.replace(microbov,method="0")
#' as(x[1:3],"data.frame")
#' 
#' ## dudi functions attempt to convert their first argument
#' ## to a data.frame; so they can be used on genind/genpop objects.
#' ## perform a PCA
#' pca1 <- dudi.pca(x, scale=FALSE, scannf=FALSE)
#' pca1
#' 
#' x <- genind2genpop(microbov,miss="chi2")
#' x <- as(x,"ktab")
#' class(x)
#' ## perform a STATIS analysis
#' statis1 <- statis(x, scannf=FALSE)
#' statis1
#' plot(statis1)
#' 
#' }
#' 
NULL





#' genind constructor
#' 
#' Constructor for \linkS4class{genind} objects.\cr The function \code{genind}
#' creates a \linkS4class{genind} object from a matrix of allelic frequency
#' where genotypes are in rows and alleles in columns. This table must have
#' correct names for rows and columns.\cr
#' 
#' The function \code{as.genind} is an alias for \code{genind} function.\cr
#' 
#' \code{is.genind} tests if an object is a valid genind object.\cr
#' 
#' Note: to get the manpage about \linkS4class{genind}, please type 'class ?
#' genind'.
#' 
#' 
#' @aliases genind-methods genind as.genind is.genind
#' @param tab A table corresponding to the @@tab slot of a genind object, with
#' individuals in rows and alleles in columns.  Its content depends on
#' \code{type} (type of marker).\cr - 'codom': table contains allele
#' frequencies (numeric values summing to 1).\cr - 'PA': table contains binary
#' values, which indicate presence(1)/absence(0) of alleles.\cr
#' @param pop a factor giving the population of each genotype in 'x'
#' @param prevcall call of an object
#' @param ploidy an integer indicating the degree of ploidy of the genotypes.
#' Beware: 2 is not an integer, but as.integer(2) is.
#' @param type a character string indicating the type of marker: 'codom' stands
#' for 'codominant' (e.g. microstallites, allozymes); 'PA' stands for
#' 'presence/absence' (e.g. AFLP).
#' @param x an object
#' @return For \code{genind} and \code{as.genind}, a genind object. For
#' \code{is.genind}, a logical.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\linkS4class{genind} class}, and \code{\link{import2genind}}
#' for importing from various types of file.\cr
#' 
#' Related classes:\cr - \linkS4class{genpop} for storing data per
#' populations\cr
#' 
#' - \linkS4class{genlight} for an efficient storage of binary SNPs
#' genotypes\cr
#' @keywords manip
#' @examples
#' 
#' data(nancycats)
#' nancycats@@loc.names
#' 
#' # isolate one marker, fca23
#' obj <- seploc(nancycats)$"fca23"
#' obj
#' 
NULL





#' Conversion to class "genlight"
#' 
#' The class \code{genlight} is a formal (S4) class for storing a genotypes of
#' binary SNPs in a compact way, using a bit-level coding scheme. New instances
#' of this class are best created using \code{new}; see the manpage of
#' \linkS4class{genlight} for more information on this point.
#' 
#' As a shortcut, conversion methods can be used to convert various objects
#' into a \linkS4class{genlight} object. Conversions can be achieved using
#' S3-style (\code{as.genlight(x)}) or S4-style (\code{as(x,"genlight"})
#' procedures. All of them call upon the constructor (\code{new}) of
#' \linkS4class{genlight} objects.
#' 
#' Conversion is currently available from the following objects: - matrix of
#' type integer/numeric - data.frame with integer/numeric data - list of
#' vectors of integer/numeric type
#' 
#' 
#' @aliases as,genlight,matrix-method as,genlight,data.frame-method
#' as,genlight,list-method as.genlight as.genlight,matrix-method
#' as.genlight,data.frame-method as.genlight,list-method
#' coerce,genlight,matrix-method coerce,genlight,data.frame-method
#' coerce,genlight,list-method
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#' @seealso Related class:\cr - \code{\linkS4class{SNPbin}}, for storing
#' individual genotypes of binary SNPs\cr
#' 
#' - \code{\linkS4class{genind}}
#' @keywords classes
#' @examples
#' 
#' \dontrun{
#' ## data to be converted
#' dat <- list(toto=c(1,1,0,0,2,2,1,2,NA), titi=c(NA,1,1,0,1,1,1,0,0), tata=c(NA,0,3, NA,1,1,1,0,0))
#' 
#' ## using the constructor
#' x1 <- new("genlight", dat)
#' x1
#' 
#' ## using 'as' methods
#' x2 <- as.genlight(dat)
#' x3 <- as(dat, "genlight")
#' 
#' identical(x1,x2)
#' identical(x1,x3)
#' }
#' 
#' 
NULL





#' genpop constructor
#' 
#' Constructor for \linkS4class{genpop} objects.\cr The function \code{genpop}
#' creates a \linkS4class{genpop} object from a matrix of alleles counts where
#' genotypes are in rows and alleles in columns. This table must have correct
#' names for rows and columns.\cr
#' 
#' The function \code{as.genpop} is an alias for \code{genpop} function.\cr
#' 
#' \code{is.genpop} tests if an object is a valid genpop object.\cr
#' 
#' Note: to get the manpage about \linkS4class{genpop}, please type 'class ?
#' genpop'.
#' 
#' 
#' @aliases genpop-methods genpop as.genpop is.genpop
#' @param tab a pop x alleles matrix which terms are numbers of alleles, i.e.
#' like in a genpop object
#' @param prevcall call of an object
#' @param ploidy an integer indicating the degree of ploidy of the genotypes.
#' Beware: 2 is not an integer, but as.integer(2) is.
#' @param type a character string indicating the type of marker: 'codom' stands
#' for 'codominant' (e.g. microstallites, allozymes); 'PA' stands for
#' 'presence/absence' (e.g. AFLP, RAPD).
#' @param x an object
#' @return For \code{genpop} and \code{as.genpop}, a genpop object. For
#' \code{is.genpop}, a logical.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\linkS4class{genpop} class}, and \code{\link{genind2genpop}}
#' for conversion from a genind to a genpop object.
#' @keywords manip
#' @examples
#' 
#' data(nancycats)
#' obj <- genind2genpop(nancycats)
#' 
#' # isolate one locus, fca77
#' obj <- seploc(obj)$"fca77"
#' obj
#' 
NULL





#' Conversion to class "SNPbin"
#' 
#' The class \linkS4class{SNPbin} is a formal (S4) class for storing a genotype
#' of binary SNPs in a compact way, using a bit-level coding scheme. New
#' instances of this class are best created using \code{new}; see the manpage
#' of \linkS4class{SNPbin} for more information on this point.
#' 
#' As a shortcut, conversion methods can be used to convert various objects
#' into a \linkS4class{SNPbin} object. Conversions can be achieved using
#' S3-style (\code{as.SNPbin(x)}) or S4-style (\code{as(x,"SNPbin"})
#' procedures. All of them call upon the constructor (\code{new}) of
#' \linkS4class{SNPbin} objects.
#' 
#' Conversion is currently available from the following objects: - integer
#' vectors - numeric vectors
#' 
#' 
#' @aliases as,SNPbin,integer-method as,SNPbin,numeric-method as.SNPbin
#' as.SNPbin,integer-method as.SNPbin,numeric-method
#' coerce,integer,SNPbin-method coerce,numeric,SNPbin-method
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#' @seealso Related class:\cr - \code{\linkS4class{SNPbin}} -
#' \code{\linkS4class{genlight}}, for storing multiple binary SNP genotypes.
#' \cr
#' @keywords classes
#' @examples
#' 
#' \dontrun{
#' ## data to be converted
#' dat <- c(1,0,0,2,1,1,1,2,2,1,1,0,0,1)
#' 
#' ## using the constructor
#' x1 <- new("SNPbin", dat)
#' x1
#' 
#' ## using 'as' methods
#' x2 <- as.SNPbin(dat)
#' x3 <- as(dat, "SNPbin")
#' 
#' identical(x1,x2)
#' identical(x1,x3)
#' }
#' 
NULL





#' Compute and optimize a-score for Discriminant Analysis of Principal
#' Components (DAPC)
#' 
#' These functions are under development. Please email the author before using
#' them for published results.
#' 
#' The Discriminant Analysis of Principal Components seeks a reduced space
#' inside which observations are best discriminated into pre-defined groups.
#' One way to assess the quality of the discrimination is looking at
#' re-assignment of individuals to their prior group, successful re-assignment
#' being a sign of strong discrimination.
#' 
#' However, when the original space is very large, ad hoc solutions can be
#' found, which discriminate very well the sampled individuals but would
#' perform poorly on new samples. In such a case, DAPC re-assignment would be
#' high even for randomly chosen clusters.  The a-score measures this bias. It
#' is computed as (Pt-Pr), where Pt is the reassignment probability using the
#' true cluster, and Pr is the reassignment probability for randomly permuted
#' clusters. A a-score close to one is a sign that the DAPC solution is both
#' strongly discriminating and stable, while low values (toward 0 or lower)
#' indicate either weak discrimination or instability of the results.
#' 
#' The a-score can serve as a criterion for choosing the optimal number of PCs
#' in the PCA step of DAPC, i.e. the number of PC maximizing the a-score. Two
#' procedures are implemented in \code{optim.a.score}. The smart procedure
#' selects evenly distributed number of PCs in a pre-defined range, compute the
#' a-score for each, and then interpolate the results using splines, predicting
#' an approximate optimal number of PCs. The other procedure (when \code{smart}
#' is FALSE) performs the computations for all number of PCs request by the
#' user. The 'optimal' number is then the one giving the highest mean a-score
#' (computed over the groups).
#' 
#' @aliases a.score optim.a.score
#' @param x a \code{dapc} object.
#' @param n.pca a vector of \code{integers} indicating the number of axes
#' retained in the Principal Component Analysis (PCA) steps of DAPC.
#' \code{nsim} DAPC will be run for each value in \code{n.pca}, unless the
#' smart approach is used (see details).
#' @param smart a \code{logical} indicating whether a smart, less
#' computer-intensive approach should be used (TRUE, default) or not (FALSE).
#' See details section.
#' @param n an \code{integer} indicating the numbers of values spanning the
#' range of \code{n.pca} to be used in the smart approach.
#' @param plot a \code{logical} indicating whether the results should be
#' displayed graphically (TRUE, default) or not (FALSE).
#' @param n.sim an \code{integer} indicating the number of simulations to be
#' performed for each number of retained PC.
#' @param n.da an \code{integer} indicating the number of axes retained in the
#' Discriminant Analysis step.
#' @param list() further arguments passed to other methods; currently unused..
#' @return === a.score ===\cr \code{a.score} returns a list with the following
#' components:\cr \item{tab}{a matrix of a-scores with groups in columns and
#' simulations in row.} \item{pop.score}{a vector giving the mean a-score for
#' each population.} \item{mean}{the overall mean a-score.}\cr
#' 
#' === optim.a.score ===\cr \code{optima.score} returns a list with the
#' following components:\cr \item{pop.score}{a list giving the mean a-score of
#' the populations for each number of retained PC (each element of the list
#' corresponds to a number of retained PCs).} \item{mean}{a vector giving the
#' overall mean a-score for each number of retained PCs.} \item{pred}{(only
#' when \code{smart} is TRUE) the predictions of the spline, given in x and y
#' coordinates.} \item{best}{the optimal number of PCs to be retained.}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso - \code{\link{find.clusters}}: to identify clusters without prior.
#' 
#' - \code{\link{dapc}}: the Discriminant Analysis of Principal Components
#' (DAPC)
#' @references Jombart T, Devillard S and Balloux F (2010) Discriminant
#' analysis of principal components: a new method for the analysis of
#' genetically structured populations. BMC Genetics11:94.
#' doi:10.1186/1471-2156-11-94
#' @keywords multivariate
NULL





#' Auxiliary functions for adegenet
#' 
#' adegenet implements a number of auxiliary procedures that might be of
#' interest for users. These include graphical tools to translate variables
#' (numeric or factors) onto a color scale, adding transparency to existing
#' colors, pre-defined color palettes, extra functions to access documentation,
#' and low-level treatment of character vectors.
#' 
#' These functions are mostly auxiliary procedures used internally in
#' adegenet.\cr
#' 
#' These items include: \itemize{ \item \code{num2col}: translates a numeric
#' vector into colors.  \item \code{fac2col}: translates a factor into colors.
#' \item \code{any2col}: translates a vector of type numeric, character or
#' factor into colors.  \item \code{transp}: adds transparency to a vector of
#' colors. Note that transparent colors are not supported on some graphical
#' devices.  \item \code{corner}: adds text to a corner of a figure.  \item
#' \code{checkType}: checks the type of markers being used in a function and
#' issues an error if appropriate.  \item \code{.rmspaces}: remove peripheric
#' spaces in a character string.  \item \code{.genlab}: generate labels in a
#' correct alphanumeric ordering.  \item \code{.readExt}: read the extension of
#' a given file.  }
#' 
#' Color palettes include: \itemize{ \item \code{bluepal}: white -> dark blue
#' \item \code{redpal}: white -> dark red \item \code{greenpal}: white -> dark
#' green \item \code{greypal}: white -> dark grey \item \code{flame}: gold ->
#' red \item \code{azur}: gold -> blue \item \code{seasun}: blue -> gold -> red
#' \item \code{lightseasun}: blue -> gold -> red (light variant) \item
#' \code{deepseasun}: blue -> gold -> red (deep variant) \item \code{spectral}:
#' red -> yellow -> blue (RColorBrewer variant) \item \code{wasp}: gold ->
#' brown -> black \item \code{funky}: many colors }
#' 
#' 
#' @name Auxiliary functions
#' @aliases checkType .rmspaces .genlab .readExt corner num2col fac2col any2col
#' transp bluepal redpal greenpal greypal flame azur seasun lightseasun
#' deepseasun spectral wasp funky
#' @docType methods
#' @param base a character string forming the base of the labels
#' @param n the number of labels to generate
#' @param text a character string to be added to the plot
#' @param posi a character matching any combinations of "top/bottom" and
#' "left/right".
#' @param inset a vector of two numeric values (recycled if needed) indicating
#' the inset, as a fraction of the plotting region.
#' @param \dots further arguments to be passed to \code{\link{text}}
#' @param x a numeric vector (for \code{num2col}) or a vector converted to a
#' factor (for \code{fac2col}).
#' @param col.pal a function generating colors according to a given palette.
#' @param reverse a logical stating whether the palette should be inverted
#' (TRUE), or not (FALSE, default).
#' @param x.min the minimal value from which to start the color scale
#' @param x.max the maximal value from which to start the color scale
#' @param na.col the color to be used for missing values (NAs)
#' @param seed a seed for R's random number generated, used to fix the random
#' permutation of colors in the palette used; if NULL, no randomization is used
#' and the colors are taken from the palette according to the ordering of the
#' levels.
#' @param col a vector of colors
#' @param alpha a numeric value between 0 and 1 representing the alpha
#' coefficient; 0: total transparency; 1: no transparency.
#' @return For \code{.genlab}, a character vector of size "n".  \code{num2col}
#' and \code{fac2col} return a vector of colors. \code{any2col} returns a list
#' with the following components: \code{$col} (a vector of colors),
#' \code{$leg.col} (colors for the legend), and \code{$leg.txt} (text for the
#' legend).
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso The R package RColorBrewer, proposing a nice selection of color
#' palettes.
#' @keywords manip
#' @examples
#' 
#' 
#' .genlab("Locus-",11)
#' 
#' ## transparent colors using "transp"
#' plot(rnorm(1000), rnorm(1000), col=transp("blue",.3), pch=20, cex=4)
#' 
#' 
#' ## numeric values to color using num2col
#' plot(1:100, col=num2col(1:100), pch=20, cex=4)
#' plot(1:100, col=num2col(1:100, col.pal=bluepal), pch=20, cex=4)
#' plot(1:100, col=num2col(1:100, col.pal=flame), pch=20, cex=4)
#' plot(1:100, col=num2col(1:100, col.pal=wasp), pch=20, cex=4)
#' plot(1:100, col=num2col(1:100, col.pal=azur,rev=TRUE), pch=20, cex=4)
#' plot(1:100, col=num2col(1:100, col.pal=spectral), pch=20, cex=4)
#' 
#' ## factor as colors using fac2col
#' dat <- cbind(c(rnorm(50,8), rnorm(100), rnorm(150,3),
#' rnorm(50,10)),c(rnorm(50,1),rnorm(100),rnorm(150,3), rnorm(50,5)))
#' fac <- rep(letters[1:4], c(50,100,150,50))
#' plot(dat, col=fac2col(fac), pch=19, cex=4)
#' plot(dat, col=transp(fac2col(fac)), pch=19, cex=4)
#' plot(dat, col=transp(fac2col(fac,seed=2)), pch=19, cex=4)
#' 
#' ## use of any2col
#' x <- factor(1:10)
#' col.info <- any2col(x, col.pal=funky)
#' plot(x, col=col.info$col, main="Use of any2col on a factor")
#' legend("bottomleft", fill=col.info$leg.col, legend=col.info$leg.txt, bg="white")
#' 
#' x <- 100:1
#' col.info <- any2col(x, col.pal=wasp)
#' barplot(x, col=col.info$col, main="Use of any2col on a numeric")
#' legend("bottomleft", fill=col.info$leg.col, legend=col.info$leg.txt, bg="white")
#' 
#' 
NULL





#' Graphics for Discriminant Analysis of Principal Components (DAPC)
#' 
#' These functions provide graphic outputs for Discriminant Analysis of
#' Principal Components (DAPC, Jombart et al. 2010). See \code{?dapc} for
#' details about this method. DAPC graphics are detailed in the DAPC tutorial
#' accessible using \code{vignette("adegenet-dapc")}.
#' 
#' These functions all require an object of class \code{dapc} (the ".dapc" can
#' be ommitted when calling the functions):\cr - \code{scatter.dapc}: produces
#' scatterplots of principal components (or 'discriminant functions'), with a
#' screeplot of eigenvalues as inset.\cr - \code{assignplot}: plot showing the
#' probabilities of assignment of individuals to the different clusters.\cr -
#' \code{compoplot}: barplot showing the probabilities of assignment of
#' individuals to the different clusters.\cr
#' 
#' See the documentation of \code{\link{dapc}} for more information about the
#' method.
#' 
#' @aliases scatter.dapc assignplot compoplot
#' @param x a \code{dapc} object.
#' @param xax,yax \code{integers} specifying which principal components of DAPC
#' should be shown in x and y axes.
#' @param grp a factor defining group membership for the individuals. The
#' scatterplot is optimal only for the default group, i.e. the one used in the
#' DAPC analysis.
#' @param col a suitable color to be used for groups. The specified vector
#' should match the number of groups, not the number of individuals.
#' @param pch a \code{numeric} indicating the type of point to be used to
#' indicate the prior group of individuals (see \code{\link{points}}
#' documentation for more details); one value is expected for each group;
#' recycled if necessary.
#' @param bg the color used for the background of the scatterplot.
#' @param solid a value between 0 and 1 indicating the alpha level for the
#' colors of the plot; 0=full transparency, 1=solid colours.
#' @param scree.da a logical indicating whether a screeplot of Discriminant
#' Analysis eigenvalues should be displayed in inset (TRUE) or not (FALSE).
#' @param scree.pca a logical indicating whether a screeplot of Principal
#' Component Analysis eigenvalues should be displayed in inset (TRUE) or not
#' (FALSE); retained axes are displayed in black.
#' @param posi.da the position of the inset of DA eigenvalues; can match any
#' combination of "top/bottom" and "left/right".
#' @param posi.pca the position of the inset of PCA eigenvalues; can match any
#' combination of "top/bottom" and "left/right".
#' @param bg.inset the color to be used as background for the inset plots.
#' @param ratio.da the size of the inset of DA eigenvalues as a proportion of
#' the current plotting region.
#' @param ratio.pca the size of the inset of PCA eigenvalues as a proportion of
#' the current plotting region.
#' @param inset.da a vector with two numeric values (recycled if needed)
#' indicating the inset to be used for the screeplot of DA eigenvalues as a
#' proportion of the current plotting region; see \code{?add.scatter} for more
#' details.
#' @param inset.pca a vector with two numeric values (recycled if needed)
#' indicating the inset to be used for the screeplot of PCA eigenvalues as a
#' proportion of the current plotting region; see \code{?add.scatter} for more
#' details.
#' @param inset.solid a value between 0 and 1 indicating the alpha level for
#' the colors of the inset plots; 0=full transparency, 1=solid colours.
#' @param onedim.filled a logical indicating whether curves should be filled
#' when plotting a single discriminant function (TRUE), or not (FALSE).
#' @param mstree a logical indicating whether a minimum spanning tree linking
#' the groups and based on the squared distances between the groups inside the
#' entire space should added to the plot (TRUE), or not (FALSE).
#' @param lwd,lty,segcol the line width, line type, and segment colour to be
#' used for the minimum spanning tree.
#' @param legend a logical indicating whether a legend for group colours should
#' added to the plot (TRUE), or not (FALSE).
#' @param posi.leg the position of the legend for group colours; can match any
#' combination of "top/bottom" and "left/right", or a set of x/y coordinates
#' stored as a list (\code{locator} can be used).
#' @param cleg a size factor used for the legend.
#' @param
#' cstar,cellipse,axesell,label,clabel,xlim,ylim,grid,addaxes,origin,include.origin,sub,csub,possub,cgrid,pixmap,contour,area
#' arguments passed to \code{\link[ade4]{s.class}}; see \code{?s.class} for
#' more informations
#' @param only.grp a \code{character} vector indicating which groups should be
#' displayed. Values should match values of \code{x$grp}. If \code{NULL}, all
#' results are displayed
#' @param subset \code{integer} or \code{logical} vector indicating which
#' individuals should be displayed. If \code{NULL}, all results are displayed
#' @param new.pred an optional list, as returned by the \code{predict} method
#' for \code{dapc} objects; if provided, the individuals with unknown groups
#' are added at the bottom of the plot. To visualize these individuals only,
#' specify \code{only.grp="unknown"}.
#' @param cex.lab a \code{numeric} indicating the size of labels.
#' @param lab a vector of characters (recycled if necessary) of labels for the
#' individuals; if left to NULL, the row names of \code{x$tab} are used.
#' @param txt.leg a character vector indicating the text to be used in the
#' legend; if not provided, group names stored in \code{x$grp} are used.
#' @param ncol an integer indicating the number of columns of the legend,
#' defaulting to 4.
#' @param posi a characther string indicating the position of the legend; can
#' match any combination of "top/bottom" and "left/right". See \code{?legend}.
#' @param list() further arguments to be passed to other functions. For
#' \code{scatter}, arguments passed to \code{points}; for \code{compoplot},
#' arguments passed to \code{barplot}.
#' @return All functions return the matched call.\cr
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso - \code{\link{dapc}}: implements the DAPC.
#' 
#' - \code{\link{find.clusters}}: to identify clusters without prior.
#' 
#' - \code{\link{dapcIllus}}: a set of simulated data illustrating the DAPC
#' 
#' - \code{\link{eHGDP}}, \code{\link{H3N2}}: empirical datasets illustrating
#' DAPC
#' @references Jombart T, Devillard S and Balloux F (2010) Discriminant
#' analysis of principal components: a new method for the analysis of
#' genetically structured populations. BMC Genetics11:94.
#' doi:10.1186/1471-2156-11-94
#' @keywords multivariate
#' @examples
#' 
#' \dontrun{
#' data(H3N2)
#' dapc1 <- dapc(H3N2, pop=H3N2$other$epid, n.pca=30,n.da=6)
#' 
#' ## defautl plot ##
#' scatter(dapc1)
#' 
#' ## showing different scatter options ##
#' ## remove internal segments and ellipses, different pch, add MStree
#' scatter(dapc1, pch=18:23, cstar=0, mstree=TRUE, lwd=2, lty=2, posi.da="topleft")
#' 
#' ## only ellipse, custom labels, use insets
#' scatter(dapc1, cell=2, pch="", cstar=0, posi.pca="topleft", posi.da="topleft", scree.pca=TRUE,
#' inset.pca=c(.01,.3), lab=paste("year\n",2001:2006), axesel=FALSE, col=terrain.colors(10))
#' 
#' ## without ellipses, use legend for groups
#' scatter(dapc1, cell=0, cstar=0, scree.da=FALSE, clab=0, cex=3,
#' solid=.4, bg="white", leg=TRUE, posi.leg="topleft")
#' 
#' ## only one axis
#' scatter(dapc1,1,1,scree.da=FALSE, legend=TRUE, solid=.4,bg="white")
#' 
#' 
#' 
#' ## example using genlight objects ##
#' ## simulate data
#' x <- glSim(50,4e3-50, 50, ploidy=2)
#' x
#' plot(x)
#' 
#' ## perform DAPC
#' dapc2 <- dapc(x, n.pca=10, n.da=1)
#' dapc2
#' 
#' ## plot results
#' scatter(dapc2, scree.da=FALSE, leg=TRUE, txt.leg=paste("group",
#' c('A','B')), col=c("red","blue"))
#' 
#' ## SNP contributions
#' loadingplot(dapc2$var.contr)
#' loadingplot(tail(dapc2$var.contr, 100), main="Loading plot - last 100 SNPs")
#' 
#' 
#' 
#' ## assignplot / compoplot ##
#' assignplot(dapc1, only.grp=2006)
#' 
#' data(microbov)
#' dapc3 <- dapc(microbov, n.pca=20, n.da=15)
#' compoplot(dapc3, lab="")
#' }
#' 
NULL





#' Simulated data illustrating the DAPC
#' 
#' Datasets illustrating the Discriminant Analysis of Principal Components
#' (DAPC, Jombart et al. submitted).\cr
#' 
#' These data were simulated using various models using Easypop (2.0.1).  The
#' \code{dapcIllus} is a list containing the following \linkS4class{genind}
#' objects:\cr - "a": island model with 6 populations \cr - "b": hierarchical
#' island model with 6 populations (3,2,1) \cr - "c": one-dimensional stepping
#' stone with 2x6 populations, and a boundary between the two sets of 6
#' populations\cr - "d": one-dimensional stepping stone with 24 populations\cr
#' 
#' See "source" for a reference providing simulation details.
#' 
#' 
#' @name dapcIllus
#' @docType data
#' @format \code{dapcIllus} is list of 4 components being all genind objects.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso - \code{\link{dapc}}: implements the DAPC.
#' 
#' - \code{\link{eHGDP}}: dataset illustrating the DAPC and
#' \code{find.clusters}.
#' 
#' - \code{\link{H3N2}}: dataset illustrating the DAPC.
#' 
#' - \code{\link{find.clusters}}: to identify clusters without prior.
#' @references Jombart, T., Devillard, S. and Balloux, F.  Discriminant
#' analysis of principal components: a new method for the analysis of
#' genetically structured populations. Submitted to \emph{Genetics}.
#' @source Jombart, T., Devillard, S. and Balloux, F.  Discriminant analysis of
#' principal components: a new method for the analysis of genetically
#' structured populations. Submitted to \emph{BMC genetics}.
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' 
#' data(dapcIllus)
#' attach(dapcIllus)
#' a # this is a genind object, like b, c, and d.
#' 
#' 
#' ## FINS CLUSTERS EX NIHILO
#' clust.a <- find.clusters(a, n.pca=100, n.clust=6)
#' clust.b <- find.clusters(b, n.pca=100, n.clust=6)
#' clust.c <- find.clusters(c, n.pca=100, n.clust=12)
#' clust.d <- find.clusters(d, n.pca=100, n.clust=24)
#' 
#' ## examin outputs
#' names(clust.a)
#' lapply(clust.a, head)
#' 
#' 
#' ## PERFORM DAPCs
#' dapc.a <- dapc(a, pop=clust.a$grp, n.pca=100, n.da=5)
#' dapc.b <- dapc(b, pop=clust.b$grp, n.pca=100, n.da=5)
#' dapc.c <- dapc(c, pop=clust.c$grp, n.pca=100, n.da=11)
#' dapc.d <- dapc(d, pop=clust.d$grp, n.pca=100, n.da=23)
#' 
#' 
#' ## LOOK AT ONE RESULT
#' dapc.a
#' summary(dapc.a)
#' 
#' ## FORM A LIST OF RESULTS FOR THE 4 DATASETS
#' lres <- list(dapc.a, dapc.b, dapc.c, dapc.d)
#' 
#' 
#' ## DRAW 4 SCATTERPLOTS
#' par(mfrow=c(2,2))
#' lapply(lres, scatter)
#' 
#' 
#' # detach data
#' detach(dapcIllus)
#' }
#' 
NULL





#' Extended HGDP-CEPH dataset
#' 
#' This dataset consists of 1350 individuals from native Human populations
#' distributed worldwide typed at 678 microsatellite loci. The original
#' HGDP-CEPH panel [1-3] has been extended by several native American
#' populations [4]. This dataset was used to illustrate the Discriminant
#' Analysis of Principal Components (DAPC, [5]).
#' 
#' 
#' @name eHGDP
#' @docType data
#' @format \code{eHGDP} is a genind object with a data frame named
#' \code{popInfo} as supplementary component (\code{eHGDP@@other$popInfo}),
#' which contains the following variables: \describe{ \item{Population: }{a
#' character vector indicating populations.} \item{Region: }{a character vector
#' indicating the geographic region of each population.} \item{Label: }{a
#' character vector indicating the correspondence with population labels used
#' in the genind object (i.e., as output by \code{pop(eHGDP)}).}
#' \item{Latitude,Longitude: }{geographic coordinates of the populations,
#' indicated as north and east degrees.} }
#' @references [1] Rosenberg NA, Pritchard JK, Weber JL, Cann HM, Kidd KK, et
#' al. (2002) Genetic structure of human populations. \emph{Science} 298:
#' 2381-2385.
#' 
#' [2] Ramachandran S, Deshpande O, Roseman CC, Rosenberg NA, Feldman MW, et
#' al. (2005) Support from the relationship of genetic and geographic distance
#' in human populations for a serial founder effect originating in Africa.
#' \emph{Proc Natl Acad Sci U S A} 102: 15942-15947.
#' 
#' [3] Cann HM, de Toma C, Cazes L, Legrand MF, Morel V, et al. (2002) A human
#' genome diversity cell line panel.  \emph{Science} 296: 261-262.
#' 
#' [4] Wang S, Lewis CM, Jakobsson M, Ramachandran S, Ray N, et al. (2007)
#' Genetic Variation and Population Structure in Native Americans. \emph{PLoS
#' Genetics} 3: e185.
#' 
#' [5] Jombart, T., Devillard, S. and Balloux, F.  Discriminant analysis of
#' principal components: a new method for the analysis of genetically
#' structured populations. Submitted to \emph{BMC genetics}.
#' @source Original panel by Human Genome Diversity Project (HGDP) and Centre
#' d'Etude du Polymorphisme Humain (CEPH). See reference [4] for Native
#' American populations.
#' 
#' This copy of the dataset was prepared by Francois Balloux.
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' ## LOAD DATA
#' data(eHGDP)
#' eHGDP
#' 
#' 
#' ## PERFORM DAPC - USE POPULATIONS AS CLUSTERS
#' ## to reproduce exactly analyses from the paper, use "n.pca=1000"
#' dapc1 <- dapc(eHGDP, all.contrib=TRUE, scale=FALSE,
#' n.pca=200, n.da=80) # takes 2 minutes
#' dapc1
#' 
#' ## (see ?dapc for details about the output)
#' 
#' 
#' 
#' ## SCREEPLOT OF EIGENVALUES
#' barplot(dapc1$eig, main="eHGDP - DAPC eigenvalues",
#' col=c("red","green","blue", rep("grey", 1000)))
#' 
#' 
#' 
#' ## SCATTERPLOTS
#' ## (!) Note: colors may be inverted with respect to [5]
#' ## as signs of principal components are arbitrary
#' ## and change from one computer to another
#' ##
#' ## axes 1-2
#' s.label(dapc1$grp.coord[,1:2], clab=0, sub="Axes 1-2")
#' par(xpd=T)
#' colorplot(dapc1$grp.coord[,1:2], dapc1$grp.coord, cex=3, add=TRUE)
#' add.scatter.eig(dapc1$eig,10,1,2, posi="bottomright", ratio=.3, csub=1.25)
#' 
#' ## axes 2-3
#' s.label(dapc1$grp.coord[,2:3], clab=0, sub="Axes 2-3")
#' par(xpd=T)
#' colorplot(dapc1$grp.coord[,2:3], dapc1$grp.coord, cex=3, add=TRUE)
#' add.scatter.eig(dapc1$eig,10,1,2, posi="bottomright", ratio=.3, csub=1.25)
#' 
#' 
#' 
#' ## MAP DAPC1 RESULTS
#' if(require(maps)){
#' 
#' xy <- cbind(eHGDP$other$popInfo$Longitude, eHGDP$other$popInfo$Latitude)
#' 
#' par(mar=rep(.1,4))
#' map(fill=TRUE, col="lightgrey")
#' colorplot(xy, -dapc1$grp.coord, cex=3, add=TRUE, trans=FALSE)
#' }
#' 
#' 
#' 
#' ## LOOK FOR OTHER CLUSTERS
#' ## to reproduce results of the reference paper, use :
#' ## grp <- find.clusters(eHGDP, max.n=50, n.pca=200, scale=FALSE)
#' ## and then
#' ## plot(grp$Kstat, type="b", col="blue")
#' 
#' grp <- find.clusters(eHGDP, max.n=30, n.pca=200,
#' scale=FALSE, n.clust=4) # takes about 2 minutes
#' names(grp)
#' 
#' ## (see ?find.clusters for details about the output)
#' 
#' 
#' 
#' ## PERFORM DAPC - USE POPULATIONS AS CLUSTERS
#' ## to reproduce exactly analyses from the paper, use "n.pca=1000"
#' dapc2 <- dapc(eHGDP, pop=grp$grp, all.contrib=TRUE,
#' scale=FALSE, n.pca=200, n.da=80) # takes around a 1 minute
#' dapc2
#' 
#' 
#' ## PRODUCE SCATTERPLOT
#' scatter(dapc2) # axes 1-2
#' scatter(dapc2,2,3) # axes 2-3
#' 
#' 
#' ## MAP DAPC2 RESULTS
#' if(require(maps)){
#' xy <- cbind(eHGDP$other$popInfo$Longitude,
#' eHGDP$other$popInfo$Latitude)
#' 
#' myCoords <- apply(dapc2$ind.coord, 2, tapply, pop(eHGDP), mean)
#' 
#' par(mar=rep(.1,4))
#' map(fill=TRUE, col="lightgrey")
#' colorplot(xy, myCoords, cex=3, add=TRUE, trans=FALSE)
#' }
#' 
#' }
#' 
NULL





#' Conversion functions from adegenet to other R packages
#' 
#' The function \code{genind2genotype} and \code{genind2hierfstat} convert a
#' \code{genind} object into, respectively, a list of genotypes (class
#' \code{genotypes}, package \code{genetics}), and a data.frame to be used by
#' the functions of the package \code{hierfstat}.
#' 
#' 
#' @aliases genind2genotype genind2hierfstat
#' @param x a \code{genind} object.
#' @param pop a factor giving the population of each individual. If NULL, it is
#' seeked in x\$pop. If NULL again, all individuals are assumed from the same
#' population.
#' @param res.type a character (if a vector, only the first element is
#' retained), indicating the type of result returned.
#' @return The function \code{genind2genotype} converts a \code{genind} object
#' into \code{genotypes} (package \code{genetics}).\cr If res.type is set to
#' "matrix" (default), the returned value is a individuals x locus matrix whose
#' columns have the class \code{genotype}. Such data can be used by
#' \code{LDheatmap} package to compute linkage disequilibrium.\cr
#' 
#' If res.type is set to "list", the returned value is a list of
#' \code{genotypes} sorted first by locus and then by population.)\cr
#' 
#' \code{genind2hierfstat} returns a data frame where individuals are in rows.
#' The first columns is a population factor (but stored as integer); each other
#' column is a locus. Genotypes are coded as integers (e.g., 44 is an
#' homozygote 4/4, 56 is an heterozygote 5/6).\cr
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{import2genind}}
#' @references Gregory Warnes and Friedrich Leisch (2007). genetics: Population
#' Genetics. R package version 1.2.1.
#' 
#' Jerome Goudet (2005). HIERFSTAT, a package for R to compute and test
#' hierarchical F-statistics. \emph{Molecular Ecology}, \bold{5}:184-186 \cr
#' 
#' Fstat (version 2.9.3). Software by Jerome Goudet.
#' http://www2.unil.ch/popgen/softwares/fstat.htm\cr
#' @keywords manip
NULL





#' F statistics for genind objects
#' 
#' \code{pairwise.fst} computes Nei's pairwise Fst between all pairs of
#' populations using a \linkS4class{genind} object. Heretozygosities are
#' weighted by group sizes (see details).
#' 
#' The function \code{fstat} is a wrapper for \code{varcomp.glob} of the
#' package \code{hierfstat}. For Fst, Fis and Fit, an alternative is offered by
#' \code{Fst} from the \code{pagas} package (see example).
#' 
#' Let \eqn{A} and \eqn{B} be two populations of population sizes \eqn{n_A} and
#' \eqn{n_B}, with expected heterozygosity (averaged over loci) \eqn{Hs(A)} and
#' \eqn{Hs(B)}, respectively. We denote \eqn{Ht} the expected heterozygosity of
#' a population pooling \eqn{A} and \eqn{B}. Then, the pairwise \eqn{Fst}
#' between \eqn{A} and \eqn{B} is computed as:\cr
#' 
#' \eqn{ Fst(A,B) = \frac{(Ht - (n_A Hs(A) + n_B Hs(B))/(n_A + n_B) )}{Ht}} \cr
#' 
#' @aliases fstat FST fst pairwise.fst
#' @param x an object of class \linkS4class{genind}.
#' @param pop a factor giving the 'population' of each individual. If NULL, pop
#' is seeked from \code{pop(x)}. Note that the term population refers in fact
#' to any grouping of individuals'.
#' @param res.type the type of result to be returned: a \code{dist} object, or
#' a symmetric matrix
#' @param truenames a logical indicating whether true labels (as opposed to
#' generic labels) should be used to name the output.
#' @param fstonly a logical stating whether only the Fst should be returned.
#' @return A vector, a matrix, or a dist object containing F statistics.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{Hs}}
#' @references Nei, M. (1973) Analysis of gene diversity in subdivided
#' populations. Proc Natl Acad Sci USA, 70: 3321-3323
#' @keywords multivariate
#' @examples
#' 
#' data(nancycats)
#' 
#' \dontrun{
#' ## pairwise Fst
#' mat.fst <- pairwise.fst(nancycats, res.type="matrix")
#' mat.fst
#' }
#' 
#' ## Fst, Fis, Fit
#' ## using hierfstat
#' if(require(hierfstat)){
#' fstat(nancycats)
#' }
#' 
#' ## using pegas
#' if(require(pegas)){
#' data(nancycats)
#' 
#' ## conversion to pegas's format
#' as.loci(nancycats)
#' 
#' ## use Fst from pegas
#' fsttab <- Fst(as.loci(nancycats))
#' 
#' ## average over loci
#' apply(fsttab, 2, mean)
#' }
#' 
NULL





#' adegenet formal class (S4) for individual genotypes
#' 
#' The S4 class \code{genind} is used to store individual genotypes.\cr It
#' contains several components described in the 'slots' section).\cr The
#' \code{summary} of a \code{genind} object invisibly returns a list of
#' component.  The function \code{.valid.genind} is for internal use.  The
#' function \code{genind} creates a genind object from a valid table of alleles
#' corresponding to the \code{@@@tab} slot.  Note that as in other S4 classes,
#' slots are accessed using @@ instead of \$.
#' 
#' 
#' @aliases genind-class print,genind-method show,genind-method
#' names,genind-method summary,genind-method .valid.genind
#' @section Slots: \describe{ \item{list("tab")}{matrix of genotypes (in rows)
#' for all alleles (in columns). The table differs depending on the
#' \code{@@type} slot:\cr - 'codom': values are frequencies ; '0' if the
#' genotype does not have the corresponding allele, '1' for an homozygote and
#' 0.5 for an heterozygte.\cr - 'PA': values are presence/absence of
#' alleles.\cr In all cases, rows and columns are given generic
#' names.}\item{:}{matrix of genotypes (in rows) for all alleles (in columns).
#' The table differs depending on the \code{@@type} slot:\cr - 'codom': values
#' are frequencies ; '0' if the genotype does not have the corresponding
#' allele, '1' for an homozygote and 0.5 for an heterozygte.\cr - 'PA': values
#' are presence/absence of alleles.\cr In all cases, rows and columns are given
#' generic names.} \item{list("loc.names")}{character vector containing the
#' real names of the loci}\item{:}{character vector containing the real names
#' of the loci} \item{list("loc.fac")}{locus factor for the columns of
#' \code{tab}}\item{:}{locus factor for the columns of \code{tab}}
#' \item{list("loc.nall")}{integer vector giving the number of alleles per
#' locus}\item{:}{integer vector giving the number of alleles per locus}
#' \item{list("all.names")}{list having one component per locus, each
#' containing a character vector of alleles names}\item{:}{list having one
#' component per locus, each containing a character vector of alleles names}
#' \item{list("call")}{the matched call}\item{:}{the matched call}
#' \item{list("ind.names")}{character vector containing the real names of the
#' individuals. Note that as Fstat does not store these names, objects
#' converted from .dat files will contain empty
#' \code{ind.names}.}\item{:}{character vector containing the real names of the
#' individuals. Note that as Fstat does not store these names, objects
#' converted from .dat files will contain empty \code{ind.names}.}
#' \item{list("ploidy")}{ an integer indicating the degree of ploidy of the
#' genotypes. Beware: 2 is not an integer, but as.integer(2) is.}\item{:}{ an
#' integer indicating the degree of ploidy of the genotypes. Beware: 2 is not
#' an integer, but as.integer(2) is.} \item{list("type")}{ a character string
#' indicating the type of marker: 'codom' stands for 'codominant' (e.g.
#' microstallites, allozymes); 'PA' stands for 'presence/absence' (e.g.
#' AFLP).}\item{:}{ a character string indicating the type of marker: 'codom'
#' stands for 'codominant' (e.g. microstallites, allozymes); 'PA' stands for
#' 'presence/absence' (e.g. AFLP).} \item{list("pop")}{(optional) factor giving
#' the population of each individual}\item{:}{(optional) factor giving the
#' population of each individual} \item{list("pop.names")}{(optional) vector
#' giving the real names of the populations}\item{:}{(optional) vector giving
#' the real names of the populations} \item{list("other")}{(optional) a list
#' containing other information}\item{:}{(optional) a list containing other
#' information} }
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{as.genind}}, \code{\link{is.genind}},
#' \code{\link{genind2genpop}}, \code{\link{genpop}},
#' \code{\link{import2genind}}, \code{\link{read.genetix}},
#' \code{\link{read.genepop}}, \code{\link{read.fstat}},
#' \code{\link{na.replace}}\cr
#' 
#' Related classes:\cr - \linkS4class{genpop} for storing data per
#' populations\cr
#' 
#' - \linkS4class{genlight} for an efficient storage of binary SNPs
#' genotypes\cr
#' @keywords classes manip multivariate
#' @examples
#' 
#' showClass("genind")
#' 
#' obj <- read.genetix(system.file("files/nancycats.gtx",package="adegenet"),missing="mean")
#' obj
#' validObject(obj)
#' summary(obj)
#' 
#' \dontrun{
#' # test inter-colonies structuration
#' if(require(hierfstat)){
#' gtest <- gstat.randtest(obj,nsim=99)
#' gtest
#' plot(gtest)
#' }
#' 
#' # perform a between-class PCA
#' pca1 <- dudi.pca(obj@@tab,scannf=FALSE,scale=FALSE)
#' pcabet1 <- between(pca1,obj@@pop,scannf=FALSE)
#' pcabet1
#' 
#' s.class(pcabet1$ls,obj@@pop,sub="Inter-class PCA",possub="topleft",csub=2)
#' add.scatter.eig(pcabet1$eig,2,xax=1,yax=2)
#' 
#' }
#' 
NULL





#' Formal class "genlight"
#' 
#' The class \code{genlight} is a formal (S4) class for storing a genotypes of
#' binary SNPs in a compact way, using a bit-level coding scheme.  This storage
#' is most efficient with haploid data, where the memory taken to represent
#' data can be reduced more than 50 times. However, \code{genlight} can be used
#' for any level of ploidy, and still remain an efficient storage mode.
#' 
#' A \code{genlight} object can be constructed from vectors of integers giving
#' the number of the second allele for each locus and each individual (see
#' 'Objects of the class genlight' below).
#' 
#' \code{genlight} stores multiple genotypes. Each genotype is stored as a
#' \linkS4class{SNPbin} object.
#' 
#' === On the subsetting using \code{[} ===
#' 
#' The function \code{[} accepts the following extra arguments: \describe{
#' \item{treatOther}{a logical stating whether elements of the \code{@@other}
#' slot should be treated as well (TRUE), or not (FALSE). If treated, elements
#' of the list are examined for a possible match of length (vectors, lists) or
#' number of rows (matrices, data frames) with the number of individuals. Those
#' who match are subsetted accordingly. Others are left as is, issuing a
#' warning unless the argument \code{quiet} is set to TRUE.} \item{quiet}{a
#' logical indicating whether warnings should be issued when trying to subset
#' components of the \code{@@other} slot which do not match the number of
#' individuals (TRUE), or not (FALSE, default). } \item{list()}{further
#' arguments passed to the genlight constructor.} }
#' 
#' @name genlight-class
#' @aliases genlight genlight-class [,genlight-method [,genlight,ANY,ANY-method
#' initialize,genlight-method show,genlight-method nLoc,genlight-method
#' nInd,genlight-method $,genlight-method $<-,genlight-method
#' names,genlight-method ploidy,genlight-method ploidy<-,genlight-method
#' locNames,genlight-method locNames<-,genlight-method indNames,genlight-method
#' indNames<-,genlight-method alleles,genlight-method alleles<-,genlight-method
#' chromosome chromosome<- chromosome,genlight-method
#' chromosome<-,genlight-method chr chr<- chr,genlight-method
#' chr<-,genlight-method position position<- position,genlight-method
#' position<-,genlight-method pop,genlight-method pop<-,genlight-method NA.posi
#' NA.posi,genlight-method other,genlight-method other<-,genlight-method
#' as.matrix.genlight as.data.frame.genlight as,matrix,genlight-method
#' as,data.frame,genlight-method as,list,genlight-method
#' coerce,matrix,genlight-method coerce,data.frame,genlight-method
#' coerce,list,genlight-method as.list.genlight cbind.genlight rbind.genlight
#' @docType class
#' @section Objects from the class genlight: \code{genlight} objects can be
#' created by calls to \code{new("genlight", ...)}, where '...' can be the
#' following arguments: \describe{ \item{list("gen")}{input genotypes, where
#' each genotype is coded as a vector of numbers of the second allele. If a
#' list, each slot of the list correspond to an individual; if a matrix or a
#' data.frame, rows correspond to individuals and columns to SNPs. If
#' individuals or loci are named in the input, these names will we stored in
#' the produced object. All individuals are expected to have the same number of
#' SNPs. Shorter genotypes are completed with NAs, issuing a warning.}
#' \item{list("ploidy")}{an optional vector of integers indicating the ploidy
#' of the genotypes. Genotypes can therefore have different ploidy. If not
#' provided, ploidy will be guessed from the data (as the maximum number of
#' second alleles in each individual).} \item{list("ind.names")}{an optional
#' vector of characters giving the labels of the genotypes.}
#' \item{list("loc.names")}{an optional vector of characters giving the labels
#' of the SNPs.} \item{list("loc.all")}{an optional vector of characters
#' indicating the alleles of each SNP; for each SNP, alleles must be coded by
#' two letters separated by '/', e.g. 'a/t' is valid, but 'a t' or 'a |t' are
#' not.} \item{list("chromosome")}{an optional factor indicating the chromosome
#' to which each SNP belongs.} \item{list("position")}{an optional vector of
#' integers indicating the position of the SNPs.} \item{list("other")}{an
#' optional list storing miscellaneous information.} }
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#' @seealso Related class:\cr - \code{\linkS4class{SNPbin}}, for storing
#' individual genotypes of binary SNPs\cr
#' 
#' - \code{\linkS4class{genind}}, for storing other types of genetic markers.
#' \cr
#' @keywords classes
#' @examples
#' 
#' \dontrun{
#' ## TOY EXAMPLE ##
#' ## create and convert data
#' dat <- list(toto=c(1,1,0,0), titi=c(NA,1,1,0), tata=c(NA,0,3, NA))
#' x <- new("genlight", dat)
#' x
#' 
#' ## examine the content of the object
#' names(x)
#' x@@gen
#' x@@gen[[1]]@@snp # bit-level coding for first individual
#' 
#' ## conversions
#' as.list(x)
#' as.matrix(x)
#' 
#' ## round trips - must return TRUE
#' identical(x, new("genlight", as.list(x))) # list
#' identical(x, new("genlight", as.matrix(x))) # matrix
#' identical(x, new("genlight", as.data.frame(x))) # data.frame
#' 
#' ## test subsetting
#' x[c(1,3)] # keep individuals 1 and 3
#' as.list(x[c(1,3)])
#' x[c(1,3), 1:2] # keep individuals 1 and 3, loci 1 and 2
#' as.list(x[c(1,3), 1:2])
#' x[c(TRUE,FALSE), c(TRUE,TRUE,FALSE,FALSE)] # same, using logicals
#' as.list(x[c(TRUE,FALSE), c(TRUE,TRUE,FALSE,FALSE)])
#' 
#' 
#' ## REAL-SIZE EXAMPLE ##
#' ## 50 genotypes of 1,000,000 SNPs
#' dat <- lapply(1:50, function(i) sample(c(0,1,NA), 1e6, prob=c(.5, .49, .01), replace=TRUE))
#' names(dat) <- paste("indiv", 1:length(dat))
#' print(object.size(dat), unit="aut") # size of the original data
#' 
#' x <- new("genlight", dat) # conversion
#' x
#' print(object.size(x), unit="au") # size of the genlight object
#' object.size(dat)/object.size(x) # conversion efficiency
#' 
#' 
#' 
#' #### cbind, rbind ####
#' a <- new("genlight", list(toto=rep(1,10), tata=rep(c(0,1), each=5), titi=c(NA, rep(1,9)) ))
#' 
#' ara <- rbind(a,a)
#' ara
#' as.matrix(ara)
#' 
#' aca <- cbind(a,a)
#' aca
#' as.matrix(aca)
#' 
#' 
#' #### subsetting @@other ####
#' x <- new("genlight", list(a=1,b=0,c=1), other=list(1:3, letters,data.frame(2:4)))
#' x
#' other(x)
#' x[2:3]
#' other(x[2:3])
#' other(x[2:3, treatOther=FALSE])
#' 
#' 
#' #### seppop ####
#' pop(x) # no population info
#' pop(x) <- c("pop1","pop1", "pop2") # set population memberships
#' pop(x)
#' seppop(x)
#' }
#' 
NULL





#' adegenet formal class (S4) for allele counts in populations
#' 
#' An object of class \code{genpop} contain alleles counts for several loci.\cr
#' It contains several components (see 'slots' section).\cr Such object is
#' obtained using \code{genind2genpop} which converts individuals genotypes of
#' known population into a \code{genpop} object.  Note that the function
#' \code{summary} of a \code{genpop} object returns a list of components.  Note
#' that as in other S4 classes, slots are accessed using @@ instead of \$.
#' 
#' 
#' @aliases genpop-class dist,genpop,ANY,ANY,ANY,missing-method
#' names,genpop-method show,genpop-method summary,genpop-method
#' @section Slots: \describe{ \item{list("tab")}{matrix of alleles counts for
#' each combinaison of population -in rows- and alleles -in columns-. Rows and
#' columns are given generic names.}\item{:}{matrix of alleles counts for each
#' combinaison of population -in rows- and alleles -in columns-. Rows and
#' columns are given generic names.} \item{list("loc.names")}{character vector
#' containing the real names of the loci}\item{:}{character vector containing
#' the real names of the loci} \item{list("loc.fac")}{locus factor for the
#' columns of \code{tab}}\item{:}{locus factor for the columns of \code{tab}}
#' \item{list("loc.nall")}{integer vector giving the number of alleles per
#' locus}\item{:}{integer vector giving the number of alleles per locus}
#' \item{list("all.names")}{list having one component per locus, each
#' containing a character vector of alleles names}\item{:}{list having one
#' component per locus, each containing a character vector of alleles names}
#' \item{list("call")}{the matched call}\item{:}{the matched call}
#' \item{list("pop.names")}{character vector containing the real names of the
#' populations}\item{:}{character vector containing the real names of the
#' populations} \item{list("ploidy")}{ an integer indicating the degree of
#' ploidy of the genotypes. Beware: 2 is not an integer, but as.integer(2)
#' is.}\item{:}{ an integer indicating the degree of ploidy of the genotypes.
#' Beware: 2 is not an integer, but as.integer(2) is.} \item{list("type")}{ a
#' character string indicating the type of marker: 'codom' stands for
#' 'codominant' (e.g. microstallites, allozymes); 'PA' stands for
#' 'presence/absence' (e.g. AFLP).}\item{:}{ a character string indicating the
#' type of marker: 'codom' stands for 'codominant' (e.g. microstallites,
#' allozymes); 'PA' stands for 'presence/absence' (e.g. AFLP).}
#' \item{list("other")}{(optional) a list containing other
#' information}\item{:}{(optional) a list containing other information} }
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{as.genpop}},
#' \code{\link{is.genpop}},\code{\link{makefreq}}, \code{\link{genind}},
#' \code{\link{import2genind}}, \code{\link{read.genetix}},
#' \code{\link{read.genepop}}, \code{\link{read.fstat}},
#' \code{\link{na.replace}}
#' @keywords classes manip multivariate
#' @examples
#' 
#' obj1 <- import2genind(system.file("files/nancycats.gen",
#' package="adegenet"))
#' obj1
#' 
#' 
#' obj2 <- genind2genpop(obj1)
#' obj2
#' 
#' \dontrun{
#' data(microsatt)
#' # use as.genpop to convert convenient count tab to genpop
#' obj3 <- as.genpop(microsatt$tab)
#' obj3
#' 
#' all(obj3@@tab==microsatt$tab)
#' all(obj3@@pop.names==rownames(microsatt$tab))
#' # it worked
#' 
#' # perform a correspondance analysis
#' obj4 <- genind2genpop(obj1,missing="chi2")
#' ca1 <- dudi.coa(as.data.frame(obj4@@tab),scannf=FALSE)
#' s.label(ca1$li,sub="Correspondance Analysis",csub=2)
#' add.scatter.eig(ca1$eig,2,xax=1,yax=2,posi="top")
#' 
#' }
#' 
NULL





#' Auxiliary functions for genlight objects
#' 
#' These functions provide facilities for usual computations using
#' \linkS4class{genlight} objects. When ploidy varies across individuals, the
#' outputs of these functions depend on whether the information units are
#' individuals, or alleles within individuals (see details).
#' 
#' These functions are:
#' 
#' - \code{glSum}: computes the sum of the number of second allele in each SNP.
#' 
#' - \code{glNA}: computes the number of missing values in each SNP.
#' 
#' - \code{glMean}: computes the mean number of second allele in each SNP.
#' 
#' - \code{glVar}: computes the variance of the number of second allele in each
#' SNP.
#' 
#' - \code{glDotProd}: computes dot products between (possibly centred/scaled)
#' vectors of individuals - uses compiled C code - used by glPca.
#' 
#' === On the unit of information ===
#' 
#' In the cases where individuals can have different ploidy, computation of
#' sums, means, etc. of allelic data depends on what we consider as a unit of
#' information.
#' 
#' To estimate e.g. allele frequencies, unit of information can be considered
#' as the allele, so that a diploid genotype contains two samples, a triploid
#' individual, three samples, etc. In such a case, all computations are done
#' directly on the number of alleles. This corresponds to \code{alleleAsUnit =
#' TRUE}.
#' 
#' However, when the focus is put on studying differences/similarities between
#' individuals, the unit of information is the individual, and all genotypes
#' possess the same information no matter what their ploidy is. In this case,
#' computations are made after standardizing individual genotypes to relative
#' allele frequencies. This corresponds to \code{alleleAsUnit = FALSE}.
#' 
#' Note that when all individuals have the same ploidy, this distinction does
#' not hold any more.
#' 
#' @aliases glSum glNA glMean glVar glDotProd
#' @param x a \linkS4class{genlight} object
#' @param alleleAsUnit a logical indicating whether alleles are considered as
#' units (i.e., a diploid genotype equals two samples, a triploid, three, etc.)
#' or whether individuals are considered as units of information.
#' @param center a logical indicating whether SNPs should be centred to mean
#' zero.
#' @param scale a logical indicating whether SNPs should be scaled to unit
#' variance.
#' @param useC a logical indicating whether compiled C code should be used
#' (TRUE) or not (FALSE, default).
#' @param parallel a logical indicating whether multiple cores -if available-
#' should be used for the computations (TRUE, default), or not (FALSE);
#' requires the package \code{parallel} to be installed (see details); this
#' option cannot be used alongside useCoption.
#' @param n.cores if \code{parallel} is TRUE, the number of cores to be used in
#' the computations; if NULL, then the maximum number of cores available on the
#' computer is used.
#' @return A numeric vector containing the requested information.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso - \code{\linkS4class{genlight}}: class of object for storing
#' massive binary SNP data.
#' 
#' - \code{\link{dapc}}: Discriminant Analysis of Principal Components.
#' 
#' - \code{\link{glPca}}: PCA for \linkS4class{genlight} objects.
#' 
#' - \code{\link{glSim}}: a simple simulator for \linkS4class{genlight}
#' objects.
#' 
#' - \code{\link{glPlot}}: plotting \linkS4class{genlight} objects.
#' @keywords multivariate
#' @examples
#' 
#' \dontrun{
#' x <- new("genlight", list(c(0,0,1,1,0), c(1,1,1,0,0,1), c(2,1,1,1,1,NA)))
#' x
#' as.matrix(x)
#' ploidy(x)
#' 
#' ## compute statistics - allele as unit ##
#' glNA(x)
#' glSum(x)
#' glMean(x)
#' 
#' ## compute statistics - individual as unit ##
#' glNA(x, FALSE)
#' glSum(x, FALSE)
#' glMean(x, FALSE)
#' 
#' ## explanation: data are taken as relative frequencies
#' temp <- as.matrix(x)/ploidy(x)
#' apply(temp,2, function(e) sum(is.na(e))) # NAs
#' apply(temp,2,sum, na.rm=TRUE) # sum
#' apply(temp,2,mean, na.rm=TRUE) # mean
#' }
#' 
NULL





#' Seasonal influenza (H3N2) HA segment data
#' 
#' The dataset \code{H3N2} consists of 1903 strains of seasonal influenza
#' (H3N2) distributed worldwide, and typed at 125 SNPs located in the
#' hemagglutinin (HA) segment. It is stored as an R object with class
#' \linkS4class{genind} and can be accessed as usual using \code{data(H3N2)}
#' (see example). These data were gathered from DNA sequences available from
#' Genbank (http://www.ncbi.nlm.nih.gov/Genbank/).
#' 
#' The data file \code{usflu.fasta} is a toy dataset also gathered from
#' Genbank, consisting of the aligned sequences of 80 seasonal influenza
#' isolates (HA segment) sampled in the US, in \code{fasta} format. This file
#' is installed alongside the package; the path to this file is automatically
#' determined by R using \code{system.file} (see example in this manpage and in
#' ?fasta2genlight) as well.
#' 
#' 
#' @name H3N2
#' @aliases H3N2 usflu usflu.fasta USflu USflu.fasta
#' @docType data
#' @format \code{H3N2} is a genind object with several data frame as
#' supplementary components (\code{H3N2@@other) slort}, which contains the
#' following items: \describe{ \item{x}{a \code{data.frame} containing
#' miscellaneous annotations of the sequences.} \item{xy}{a matrix with two
#' columns indicating the geographic coordinates of the strains, as longitudes
#' and latitudes.} \item{epid}{a character vector indicating the epidemic of
#' the strains.} }
#' @references Jombart, T., Devillard, S. and Balloux, F. Discriminant analysis
#' of principal components: a new method for the analysis of genetically
#' structured populations. Submitted to \emph{BMC genetics}.
#' @source This dataset was prepared by Thibaut Jombart
#' (t.jombart@@imperia.ac.uk), from annotated sequences available on Genbank
#' (http://www.ncbi.nlm.nih.gov/Genbank/).
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' #### H3N2 ####
#' ## LOAD DATA
#' data(H3N2)
#' H3N2
#' 
#' ## set population to yearly epidemics
#' pop(H3N2) <- factor(H3N2$other$epid)
#' 
#' 
#' 
#' ## PERFORM DAPC - USE POPULATIONS AS CLUSTERS
#' ## to reproduce exactly analyses from the paper, use "n.pca=1000"
#' dapc1 <- dapc(H3N2, all.contrib=TRUE, scale=FALSE, n.pca=150, n.da=5)
#' dapc1
#' 
#' ## (see ?dapc for details about the output)
#' 
#' 
#' ## SCREEPLOT OF EIGENVALUES
#' barplot(dapc1$eig, main="H3N2 - DAPC eigenvalues")
#' 
#' 
#' ## SCATTERPLOT (axes 1-2)
#' scatter(dapc1, posi.da="topleft", cstar=FALSE, cex=2, pch=17:22,
#' solid=.5, bg="white")
#' 
#' 
#' 
#' 
#' #### usflu.fasta ####
#' myPath <- system.file("files/usflu.fasta",package="adegenet")
#' myPath
#' 
#' ## extract SNPs from alignments using fasta2genlight
#' ## see ?fasta2genlight for more details
#' obj <- fasta2genlight(myPath, chunk=10) # process 10 sequences at a time
#' obj
#' }
#' 
NULL





#' Importing data from several softwares to a genind object
#' 
#' Their are several ways to import genotype data to a \linkS4class{genind}
#' object: i) from a data.frame with a given format (see
#' \code{\link{df2genind}}), ii) from a file with a recognized extension, or
#' iii) from an alignement of sequences (see \code{\link{DNAbin2genind}}).\cr
#' 
#' The function \code{import2genind} detects the extension of the file given in
#' argument and seeks for an appropriate import function to create a
#' \code{genind} object.\cr Current recognized formats are :\cr - GENETIX files
#' (.gtx) \cr - Genepop files (.gen) \cr - Fstat files (.dat) \cr - STRUCTURE
#' files (.str or .stru) \cr
#' 
#' There are 3 treatments for missing values: \cr - NA: kept as NA.\cr
#' 
#' - 0: allelic frequencies are set to 0 on all alleles of the concerned locus.
#' Recommended for a PCA on compositionnal data.\cr
#' 
#' - "mean": missing values are replaced by the mean frequency of the
#' corresponding allele, computed on the whole set of individuals. Recommended
#' for a centred PCA.\cr
#' 
#' Beware: same data in different formats are not expected to produce exactly
#' the same \code{genind} objects.\cr For instance, conversions made by GENETIX
#' to Fstat may change the the sorting of the genotypes; GENETIX stores
#' individual names whereas Fstat does not; Genepop chooses a sample's name
#' from the name of its last genotype; etc.
#' 
#' @aliases import2genind import2genind
#' @param file a character string giving the path to the file to convert, with
#' the appropriate extension.
#' @param missing can be NA, 0 or "mean". See details section.
#' @param quiet logical stating whether a conversion message must be printed
#' (TRUE,default) or not (FALSE).
#' @param \dots other arguments passed to the appropriate 'read' function
#' (currently passed to \code{read.structure})
#' @return an object of the class \code{genind}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{import2genind}}, \code{\link{read.genetix}},
#' \code{\link{read.fstat}}, \code{\link{read.structure}},
#' \code{\link{read.genepop}}
#' @references Belkhir K., Borsa P., Chikhi L., Raufaste N. & Bonhomme F.
#' (1996-2004) GENETIX 4.05, logiciel sous Windows TM pour la génétique des
#' populations. Laboratoire Génome, Populations, Interactions, CNRS UMR 5000,
#' Université de Montpellier II, Montpellier (France). \cr
#' 
#' Pritchard, J.; Stephens, M. & Donnelly, P. (2000) Inference of population
#' structure using multilocus genotype data. \emph{Genetics}, \bold{155}:
#' 945-959
#' 
#' Raymond M. & Rousset F, (1995). GENEPOP (version 1.2): population genetics
#' software for exact tests and ecumenicism. \emph{J. Heredity},
#' \bold{86}:248-249 \cr
#' 
#' Fstat (version 2.9.3). Software by Jerome Goudet.
#' http://www2.unil.ch/popgen/softwares/fstat.htm\cr
#' 
#' Excoffier L. & Heckel G.(2006) Computer programs for population genetics
#' data analysis: a survival guide \emph{Nature}, \bold{7}: 745-758
#' @keywords manip
#' @examples
#' 
#' import2genind(system.file("files/nancycats.gtx",
#' package="adegenet"))
#' 
#' import2genind(system.file("files/nancycats.dat",
#' package="adegenet"))
#' 
#' import2genind(system.file("files/nancycats.gen",
#' package="adegenet"))
#' 
#' import2genind(system.file("files/nancycats.str",
#' package="adegenet"), onerowperind=FALSE, n.ind=237, n.loc=9, col.lab=1, col.pop=2, ask=FALSE)
#' 
NULL





#' Likelihood-based estimation of inbreeding
#' 
#' The function \code{inbreeding} estimates the inbreeding coefficient of an
#' individuals (F) by computing its likelihood function. It can return either
#' the density of probability of F, or a sample of F values from this
#' distribution. This operation is performed for all the individuals of a
#' \linkS4class{genind} object. Any ploidy greater than 1 is acceptable.
#' 
#' Let \eqn{F} denote the inbreeding coefficient, defined as the probability
#' for an individual to inherit two identical alleles from a single ancestor.
#' 
#' Let \eqn{p_i} refer to the frequency of allele \eqn{i} in the population.
#' Let \eqn{h} be an variable which equates 1 if the individual is homozygote,
#' and 0 otherwise. For one locus, the probability of being homozygote is
#' computed as:
#' 
#' \eqn{ F + (1-F) \sum_i p_i^2}
#' 
#' The probability of being heterozygote is: \eqn{1 - (F + (1-F) \sum_i p_i^2)}
#' 
#' The likelihood of a genotype is defined as the probability of being the
#' observed state (homozygote or heterozygote). In the case of multilocus
#' genotypes, log-likelihood are summed over the loci.
#' 
#' @aliases inbreeding
#' @param x an object of class \linkS4class{genind}.
#' @param pop a factor giving the 'population' of each individual. If NULL, pop
#' is seeked from \code{pop(x)}. Note that the term population refers in fact
#' to any grouping of individuals'.
#' @param truenames a logical indicating whether true names should be used
#' (TRUE, default) instead of generic labels (FALSE); used if res.type is
#' "matrix".
#' @param res.type a character string matching "sample", "function", or
#' "estimate" specifying whether the output should be a function giving the
#' density of probability of F values ("function"), the maximum likelihood
#' estimate of F from this distribution ("estimate"), or a sample of F values
#' taken from this distribution ("sample", default).
#' @param N an integer indicating the size of the sample to be taken from the
#' distribution of F values.
#' @param M an integer indicating the number of different F values to be used
#' to generate the sample. Values larger than N are recommended to avoid poor
#' sampling of the distribution.
#' @return A named list with one component for each individual, each of which
#' is a function or a vector of sampled F values (see \code{res.type}
#' argument).
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}\cr Zhian N.
#' Kamvar\cr
#' @seealso \code{\link{Hs}}: computation of expected heterozygosity.
#' @examples
#' 
#' \dontrun{
#' ## cattle breed microsatellite data
#' data(microbov)
#' 
#' ## isolate Lagunaire breed
#' lagun <- seppop(microbov)$Lagunaire
#' 
#' ## estimate inbreeding - return sample of F values
#' Fsamp <- inbreeding(lagun, N=30)
#' 
#' ## plot the first 10 results
#' invisible(sapply(Fsamp[1:10], function(e) plot(density(e), xlab="F",
#' xlim=c(0,1), main="Density of the sampled F values")))
#' 
#' ## compute means for all individuals
#' Fmean=sapply(Fsamp, mean)
#' hist(Fmean, col="orange", xlab="mean value of F",
#' main="Distribution of mean F across individuals")
#' 
#' ## estimate inbreeding - return proba density functions
#' Fdens <- inbreeding(lagun, res.type="function")
#' 
#' ## view function for the first individual
#' Fdens[[1]]
#' 
#' ## plot the first 10 functions
#' invisible(sapply(Fdens[1:10], plot, ylab="Density",
#' main="Density of probability of F values"))
#' 
#' ## estimate inbreeding - return maximum likelihood estimates
#' Fest <- inbreeding(lagun, res.type = "estimate")
#' mostInbred <- which.max(Fest)
#' plot(Fdens[[mostInbred]], ylab = "Density", xlab = "F",
#'      main = paste("Probability density of F values\nfor", names(mostInbred)))
#' abline(v = Fest[mostInbred], col = "red", lty = 2)
#' legend("topright", legend = "MLE", col = "red", lty = 2)
#' 
#' ## note that estimates and average samples are likely to be different.
#' plot(Fest, ylab = "F", col = "blue",
#'      main = "comparison of MLE and average sample estimates of F")
#' points(Fmean, pch = 2, col = "red")
#' arrows(x0 = 1:length(Fest), y0 = Fest, 
#'        y1 = Fmean, x1 = 1:length(Fest), length = 0.125)
#' legend("topleft", legend = c("estimate", "sample"), col = c("blue", "red"),
#'        pch = c(1, 2), title = "res.type")
#' }
#' 
NULL





#' Assess polymorphism in genind/genpop objects
#' 
#' The simple function \code{isPoly} can be used to check which loci are
#' polymorphic, or alternatively to check which alleles give rise to
#' polymorphism.
#' 
#' 
#' @name isPoly-methods
#' @aliases isPoly isPoly-methods isPoly,genind-method isPoly,genpop-method
#' @docType methods
#' @param x a \linkS4class{genind} and \linkS4class{genpop} object
#' @param by a character being "locus" or "allele", indicating whether results
#' should indicate polymorphic loci ("locus"), or alleles giving rise to
#' polymorphism ("allele").
#' @param thres a numeric value giving the minimum frequency of an allele
#' giving rise to polymorphism (defaults to 0.01).
#' @return A vector of logicals.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords methods manip
#' @examples
#' 
#' \dontrun{
#' data(nancycats)
#' isPoly(nancycats,by="loc", thres=0.1)
#' isPoly(nancycats[1:3],by="loc", thres=0.1)
#' genind2df(nancycats[1:3])
#' }
#' 
NULL





#' Microsatellites genotypes of 15 cattle breeds
#' 
#' This data set gives the genotypes of 704 cattle individuals for 30
#' microsatellites recommended by the FAO. The individuals are divided into two
#' countries (Afric, France), two species (Bos taurus, Bos indicus) and 15
#' breeds. Individuals were chosen in order to avoid pseudoreplication
#' according to their exact genealogy.
#' 
#' 
#' @name microbov
#' @docType data
#' @format \code{microbov} is a genind object with 3 supplementary components:
#' \describe{ \item{coun}{a factor giving the country of each individual (AF:
#' Afric; FR: France).} \item{breed}{a factor giving the breed of each
#' individual.} \item{spe}{is a factor giving the species of each individual
#' (BT: Bos taurus; BI: Bos indicus).} }
#' @references Lalo\"e D., Jombart T., Dufour A.-B. and Moazami-Goudarzi K.
#' (2007) Consensus genetic structuring and typological value of markers using
#' Multiple Co-Inertia Analysis. \emph{Genetics Selection Evolution}.
#' \bold{39}: 545--567.
#' @source Data prepared by Katayoun Moazami-Goudarzi and Denis Lalo\"e (INRA,
#' Jouy-en-Josas, France)
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' data(microbov)
#' microbov
#' summary(microbov)
#' 
#' # make Y, a genpop object
#' Y <- genind2genpop(microbov)
#' 
#' # make allelic frequency table
#' temp <- makefreq(Y,missing="mean")
#' X <- temp$tab
#' nsamp <- temp$nobs
#' 
#' # perform 1 PCA per marker 
#' 
#' kX <- ktab.data.frame(data.frame(X),Y@@loc.nall)
#' 
#' kpca <- list()
#' for(i in 1:30) {kpca[[i]] <- dudi.pca(kX[[i]],scannf=FALSE,nf=2,center=TRUE,scale=FALSE)}
#' 
#' 
#' sel <- sample(1:30,4)
#' col = rep('red',15)
#' col[c(2,10)] = 'darkred'
#' col[c(4,12,14)] = 'deepskyblue4'
#' col[c(8,15)] = 'darkblue'
#' 
#' # display %PCA
#' par(mfrow=c(2,2))
#' for(i in sel) {
#' s.multinom(kpca[[i]]$c1,kX[[i]],n.sample=nsamp[,i],coulrow=col,sub=Y@@loc.names[i])
#' add.scatter.eig(kpca[[i]]$eig,3,xax=1,yax=2,posi="top")
#' }
#' 
#' # perform a Multiple Coinertia Analysis
#' kXcent <- kX
#' for(i in 1:30) kXcent[[i]] <- as.data.frame(scalewt(kX[[i]],center=TRUE,scale=FALSE))
#' mcoa1 <- mcoa(kXcent,scannf=FALSE,nf=3, option="uniform")
#' 
#' # coordinated %PCA
#' mcoa.axes <- split(mcoa1$axis, Y@@loc.fac)
#' mcoa.coord <- split(mcoa1$Tli,mcoa1$TL[,1])
#' var.coord <- lapply(mcoa.coord,function(e) apply(e,2,var))
#' 
#' par(mfrow=c(2,2))
#' for(i in sel) {
#' s.multinom(mcoa.axes[[i]][,1:2],kX[[i]],n.sample=nsamp[,i],coulrow=col,sub=Y@@loc.names[i])
#' add.scatter.eig(var.coord[[i]],2,xax=1,yax=2,posi="top")
#' }
#' 
#' # reference typology
#' par(mfrow=c(1,1))
#' s.label(mcoa1$SynVar,lab=microbov@@pop.names,sub="Reference typology",csub=1.5)
#' add.scatter.eig(mcoa1$pseudoeig,nf=3,xax=1,yax=2,posi="top")
#' 
#' # typologial values
#' tv <- mcoa1$cov2
#' tv <- apply(tv,2,function(c) c/sum(c))*100
#' rownames(tv) <- Y@@loc.names
#' tv <- tv[order(Y@@loc.names),]
#' 
#' par(mfrow=c(3,1),mar=c(5,3,3,4),las=3)
#' for(i in 1:3){
#' barplot(round(tv[,i],3),ylim=c(0,12),yaxt="n",main=paste("Typological value -
#' structure",i))
#' axis(side=2,at=seq(0,12,by=2),labels=paste(seq(0,12,by=2),"%"),cex=3)
#' abline(h=seq(0,12,by=2),col="grey",lty=2)
#' }
#' }
#' 
NULL





#' Replace missing values (NA) from an object
#' 
#' The generic function \code{na.replace} replaces NA in an object by
#' appropriate values as defined by the argument \code{method}.\cr
#' 
#' Methods are defined for \linkS4class{genind} and \linkS4class{genpop}
#' objects.
#' 
#' The argument "method" have the following effects:\cr - "0": missing values
#' are set to "0". An entity (individual or population) that is not type on a
#' locus has zeros for all alleles of that locus.\cr
#' 
#' - "mean": missing values are set to the mean of the concerned allele,
#' computed on all available observations (without distinction of
#' population).\cr
#' 
#' - "chi2": if a population is not typed for a marker, the corresponding count
#' is set to that of a theoretical count in of a Chi-squared test. This is
#' obtained by the product of the sums of both margins divided by the total
#' number of alleles.
#' 
#' @name na.replace-methods
#' @aliases na.replace na.replace-methods na.replace,genind-method
#' na.replace,genpop-method
#' @docType methods
#' @param x a \linkS4class{genind} and \linkS4class{genpop} object
#' @param method a character string: can be "0" or "mean" for
#' \linkS4class{genind} objects, and "0" or "chi2" for \linkS4class{genpop}
#' objects.
#' @param quiet logical stating whether a message should be printed
#' (TRUE,default) or not (FALSE).
#' @return A \linkS4class{genind} and \linkS4class{genpop} object without
#' missing values.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords methods manip
#' @examples
#' 
#' \dontrun{
#' data(nancycats)
#' 
#' obj1 <- genind2genpop(nancycats)
#' # note missing data in this summary
#' summary(obj1)
#' 
#' # NA are all in pop 17 and marker fca45
#' which(is.na(obj1$tab),TRUE)
#' truenames(obj1)[17,]
#' 
#' # replace missing values
#' obj2 <- na.replace(obj1,"chi2")
#' obj2$loc.names
#' 
#' # missing values where replaced
#' truenames(obj2)[,obj2$loc.fac=="L4"]
#' }
#' 
NULL





#' Microsatellites genotypes of 237 cats from 17 colonies of Nancy (France)
#' 
#' This data set gives the genotypes of 237 cats (\emph{Felis catus} L.) for 9
#' microsatellites markers. The individuals are divided into 17 colonies whose
#' spatial coordinates are also provided.
#' 
#' 
#' @name nancycats
#' @docType data
#' @format \code{nancycats} is a genind object with spatial coordinates of the
#' colonies as a supplementary components (@@xy). Beware: these coordinates are
#' given for the true names (stored in @@pop.names) and not for the generic
#' names (used in @@pop).
#' @references Devillard, S.; Jombart, T. & Pontier, D. Disentangling spatial
#' and genetic structure of stray cat (\emph{Felis catus} L.) colonies in urban
#' habitat using: not all colonies are equal. submitted to \emph{Molecular
#' Ecology}
#' @source Dominique Pontier (UMR CNRS 5558, University Lyon1, France)
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' data(nancycats)
#' nancycats
#' 
#' # summary's results are stored in x
#' x <- summary(nancycats)
#' 
#' # some useful graphics
#' barplot(x$loc.nall,ylab="Alleles numbers",main="Alleles numbers
#' per locus")
#' 
#' plot(x$pop.eff,x$pop.nall,type="n",xlab="Sample size",ylab="Number of alleles")
#' text(x$pop.eff,y=x$pop.nall,lab=names(x$pop.nall))
#' 
#' par(las=3)
#' barplot(table(nancycats@@pop),ylab="Number of genotypes",main="Number of genotypes per colony")
#' 
#' # are cats structured among colonies ?
#' if(require(hierfstat)){
#' 
#' gtest <- gstat.randtest(nancycats,nsim=99)
#' gtest
#' plot(gtest)
#' 
#' 
#' dat <- genind2hierfstat(nancycats)
#' 
#' Fstat <- varcomp.glob(dat$pop,dat[,-1])
#' Fstat
#' }
#' }
#' 
NULL





#' Convert objects with obsolete classe into new objects
#' 
#' Adegenet classes changed from S3 to S4 types starting from version 1.1-0.
#' \code{old2new} has two methods for genind and genpop objects, so that old
#' adegenet objects can be retrieved and used in recent versions.
#' 
#' Optional content but \code{$pop} and \code{$pop.names} will not be
#' converted. These are to be coerced into a list and set in the \code{@@other}
#' slot of the new object.
#' 
#' @name old2new
#' @aliases old2new old2new,ANY-method old2new-methods old2new,genind-method
#' old2new,genpop-method
#' @docType methods
#' @param object a genind or genpop object in S3 version, i.e. prior
#' adegenet\_1.1-0
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords methods classes manip
NULL





#' Compute the proportion of typed elements
#' 
#' The generic function \code{propTyped} is devoted to investigating the
#' structure of missing data in adegenet objects.\cr
#' 
#' Methods are defined for \linkS4class{genind} and \linkS4class{genpop}
#' objects. They can return the proportion of available (i.e. non-missing) data
#' per individual/population, locus, or the combination of both in with case
#' the matrix indicates which entity (individual or population) was typed on
#' which locus.
#' 
#' When \code{by} is set to "both", the result is a matrix of binary data with
#' entities in rows (individuals or populations) and markers in columns. The
#' values of the matrix are 1 for typed data, and 0 for NA.
#' 
#' @name propTyped-methods
#' @aliases propTyped propTyped-methods propTyped,genind-method
#' propTyped,genpop-method
#' @docType methods
#' @param x a \linkS4class{genind} and \linkS4class{genpop} object
#' @param by a character being "ind","loc", or "both" for \linkS4class{genind}
#' object and "pop","loc", or "both" for \linkS4class{genpop} object. It
#' specifies whether proportion of typed data are provided by entity
#' ("ind"/"pop"), by locus ("loc") or both ("both"). See details.
#' @return A vector of proportion (when \code{by} equals "ind", "pop", or
#' "loc"), or a matrix of binary data (when \code{by} equals "both")
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords methods manip
#' @examples
#' 
#' \dontrun{
#' data(nancycats)
#' propTyped(nancycats,by="loc")
#' propTyped(genind2genpop(nancycats),by="both")
#' }
#' 
NULL





#' Microsatellites genotypes of 335 chamois (Rupicapra rupicapra) from the
#' Bauges mountains (France)
#' 
#' This data set contains the genotypes of 335 chamois (\emph{Rupicapra
#' rupicapra}) from the Bauges mountains, in France. No prior clustering about
#' individuals is known. Each genotype is georeferenced. These data also
#' contain a raster map of elevation of the sampling area.
#' 
#' 
#' @name rupica
#' @docType data
#' @format \code{rupica} is a genind object with 3 supplementary components
#' inside the \code{@@other} slot: \describe{ \item{xy}{a matrix containing the
#' spatial coordinates of the genotypes.} \item{mnt}{a raster map of elevation,
#' with the \code{asc} format from the \code{adehabitat} package.}
#' \item{showBauges}{a function to display the map of elevation with an
#' appropriate legend (use \code{showBauges()}).} }
#' @references Cassar S (2008) Organisation spatiale de la variabilité
#' génétique et phénotypique a l'échelle du paysage: le cas du chamois et du
#' chevreuil, en milieu de montagne. PhD Thesis. University Claude Bernard -
#' Lyon 1, France. \cr
#' 
#' Cassar S, Jombart T, Loison A, Pontier D, Dufour A-B, Jullien J-M, Chevrier
#' T, Maillard D. Spatial genetic structure of Alpine chamois (\emph{Rupicapra
#' rupicapra}): a consequence of landscape features and social factors?
#' submitted to \emph{Molecular Ecology}.
#' @source Daniel Maillard, 'Office National de la Chasse et de la Faune
#' Sauvage' (ONCFS), France.
#' @keywords datasets
#' @examples
#' 
#' if(require(adehabitat) && require(spdep)){
#' 
#' data(rupica)
#' rupica
#' 
#' ## see the sampling area
#' showBauges <- rupica$other$showBauges
#' showBauges()
#' points(rupica$other$xy,col="red")
#' 
#' ## perform a sPCA
#' spca1 <- spca(rupica,type=5,d1=0,d2=2300,plot=FALSE,scannf=FALSE,nfposi=2,nfnega=0)
#' barplot(spca1$eig,col=rep(c("black","grey"),c(2,100)),main="sPCA eigenvalues")
#' screeplot(spca1,main="sPCA eigenvalues: decomposition")
#' 
#' ## data visualization
#' showBauges(,labcex=1)
#' s.value(spca1$xy,spca1$ls[,1],add.p=TRUE,csize=.5)
#' add.scatter.eig(spca1$eig,1,1,1,posi="topleft",sub="Eigenvalues")
#' 
#' showBauges(,labcex=1)
#' s.value(spca1$xy,spca1$ls[,2],add.p=TRUE,csize=.5)
#' add.scatter.eig(spca1$eig,2,2,2,posi="topleft",sub="Eigenvalues")
#' 
#' rupica$other$showBauges()
#' colorplot(spca1$xy,spca1$li,cex=1.5,add.plot=TRUE)
#' 
#' \dontrun{
#' ## global and local tests
#' Gtest <- global.rtest(rupica@@tab,spca1$lw,nperm=999)
#' Gtest
#' plot(Gtest)
#' Ltest <- local.rtest(rupica@@tab,spca1$lw,nperm=999)
#' Ltest
#' plot(Ltest)
#' }
#' }
#' 
NULL





#' Compute scaled allele frequencies
#' 
#' The generic function \code{scaleGen} is an analogue to the \code{scale}
#' function, but is designed with further arguments giving scaling options.\cr
#' 
#' Methods are defined for \linkS4class{genind} and \linkS4class{genpop}
#' objects.  Both return data.frames of scaled allele frequencies.
#' 
#' 
#' @name scaleGen-methods
#' @aliases scaleGen scaleGen-methods scaleGen,genind-method
#' scaleGen,genpop-method
#' @docType methods
#' @param x a \linkS4class{genind} and \linkS4class{genpop} object
#' @param center a logical stating whether alleles frequencies should be
#' centred to mean zero (default to TRUE). Alternatively, a vector of numeric
#' values, one per allele, can be supplied: these values will be substracted
#' from the allele frequencies.
#' @param scale a logical stating whether alleles frequencies should be scaled
#' (default to TRUE). Alternatively, a vector of numeric values, one per
#' allele, can be supplied: these values will be substracted from the allele
#' frequencies.
#' @param truenames a logical indicating whether true labels (as opposed to
#' generic labels) should be used to name the output.
#' @param missing a character giving the treatment for missing values. Can be
#' "NA", "0" or "mean"
#' @return A matrix of scaled allele frequencies with genotypes
#' (\linkS4class{genind}) or populations in (\linkS4class{genpop}) in rows and
#' alleles in columns.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords methods manip
#' @examples
#' 
#' \dontrun{
#' ## load data
#' data(microbov)
#' obj <- genind2genpop(microbov)
#' 
#' ## compare different scaling
#' X1 <- scaleGen(obj)
#' X2 <- scaleGen(obj,met="bin")
#' 
#' ## compute PCAs
#' pcaObj <- dudi.pca(obj,scale=FALSE,scannf=FALSE) # pca with no scaling
#' pcaX1 <- dudi.pca(X1,scale=FALSE,scannf=FALSE,nf=100) # pca with usual scaling
#' pcaX2 <- dudi.pca(X2,scale=FALSE,scannf=FALSE,nf=100) # pca with scaling for binomial variance
#' 
#' ## get the loadings of alleles for the two scalings
#' U1 <- pcaX1$c1
#' U2 <- pcaX2$c1
#' 
#' 
#' ## find an optimal plane to compare loadings
#' ## use a procustean rotation of loadings tables
#' pro1 <- procuste(U1,U2,nf=2)
#' 
#' ## graphics
#' par(mfrow=c(2,2))
#' # eigenvalues
#' barplot(pcaObj$eig,main="Eigenvalues\n no scaling")
#' barplot(pcaX1$eig,main="Eigenvalues\n usual scaling")
#' barplot(pcaX2$eig,main="Eigenvalues\n 'binomial' scaling")
#' # differences between loadings of alleles
#' s.match(pro1$scor1,pro1$scor2,clab=0,sub="usual -> binom (procustean rotation)")
#' 
#' }
#' 
NULL





#' Select genotypes of well-represented populations
#' 
#' The function \code{selPopSize} checks the sample size of each population in
#' a \linkS4class{genind} object and keeps only genotypes of populations having
#' a given minimum size.
#' 
#' 
#' @name selPopSize
#' @aliases selPopSize selPopSize-methods selPopSize,ANY-method
#' selPopSize,genind-method
#' @docType methods
#' @param x a \linkS4class{genind} object
#' @param pop a vector of characters or a factor giving the population of each
#' genotype in 'x'. If not provided, seeked from x\$pop.
#' @param nMin the minimum sample size for a population to be retained. Samples
#' sizes strictly less than \code{nMin} will be discarded, those equal to or
#' greater than \code{nMin} are kept.
#' @return A \linkS4class{genind} object.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{seploc}}, \code{\link{repool}}
#' @keywords manip
#' @examples
#' 
#' \dontrun{
#' data(microbov)
#' 
#' table(pop(microbov))
#' obj <- selPopSize(microbov, n=50)
#' 
#' obj
#' table(pop(obj))
#' }
#' 
NULL





#' Separate data per locus
#' 
#' The function \code{seploc} splits an object (\linkS4class{genind},
#' \linkS4class{genpop} or \linkS4class{genlight}) by marker. For
#' \linkS4class{genind} and \linkS4class{genpop} objects, the method returns a
#' list of objects whose components each correspond to a marker. For
#' \linkS4class{genlight} objects, the methods returns blocks of SNPs.
#' 
#' 
#' @name seploc
#' @aliases seploc seploc-methods seploc,ANY-method seploc,genind-method
#' seploc,genpop-method seploc,genlight-method
#' @docType methods
#' @param x a \linkS4class{genind} or a \linkS4class{genpop} object.
#' @param truenames a logical indicating whether true names should be used
#' (TRUE, default) instead of generic labels (FALSE).
#' @param res.type a character indicating the type of returned results, a
#' genind or genpop object (default) or a matrix of data corresponding to the
#' 'tab' slot.
#' @param n.block an integer indicating the number of blocks of SNPs to be
#' returned.
#' @param block.size an integer indicating the size (in number of SNPs) of the
#' blocks to be returned.
#' @param random should blocks be formed of contiguous SNPs, or should they be
#' made or randomly chosen SNPs.
#' @param parallel a logical indicating whether multiple cores -if available-
#' should be used for the computations (TRUE, default), or not (FALSE);
#' requires the package \code{parallel} to be installed.
#' @param n.cores if \code{parallel} is TRUE, the number of cores to be used in
#' the computations; if NULL, then the maximum number of cores available on the
#' computer is used.
#' @return The function \code{seploc} returns an list of objects of the same
#' class as the initial object, or a list of matrices similar to x\$tab.\cr
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{seppop}}, \code{\link{repool}}
#' @keywords manip
#' @examples
#' 
#' \dontrun{
#' ## example on genind objects
#' data(microbov)
#' 
#' # separate all markers
#' obj <- seploc(microbov)
#' names(obj)
#' 
#' obj$INRA5
#' 
#' 
#' ## example on genlight objects
#' x <- glSim(100, 1000, 0, ploidy=2) # simulate data
#' x <- x[,order(glSum(x))] # reorder loci by frequency of 2nd allele
#' glPlot(x, main="All data") # plot data
#' foo <- seploc(x, n.block=3) # form 3 blocks
#' foo
#' glPlot(foo[[1]], main="1st block") # plot 1st block
#' glPlot(foo[[2]], main="2nd block") # plot 2nd block
#' glPlot(foo[[3]], main="3rd block") # plot 3rd block
#' 
#' foo <- seploc(x, block.size=600, random=TRUE) # split data, randomize loci
#' foo # note the different block sizes
#' glPlot(foo[[1]])
#' }
#' 
NULL





#' Separate genotypes per population
#' 
#' The function \code{seppop} splits a \linkS4class{genind} or a
#' \linkS4class{genlight} object by population, returning a list of objects
#' whose components each correspond to a population.\cr
#' 
#' For \linkS4class{genind} objects, the output can either be a list of
#' \linkS4class{genind} (default), or a list of matrices corresponding to the
#' \code{@@tab} slot.
#' 
#' 
#' @name seppop
#' @aliases seppop seppop-methods seppop,ANY-method seppop,genind-method
#' seppop,genlight-method
#' @docType methods
#' @param x a \linkS4class{genind} object
#' @param pop a factor giving the population of each genotype in 'x'. If not
#' provided, seeked in x\$pop.
#' @param truenames a logical indicating whether true names should be used
#' (TRUE, default) instead of generic labels (FALSE); used if res.type is
#' "matrix".
#' @param res.type a character indicating the type of returned results, a list
#' of \linkS4class{genind} object (default) or a matrix of data corresponding
#' to the 'tab' slots.
#' @param drop a logical stating whether alleles that are no longer present in
#' a subset of data should be discarded (TRUE) or kept anyway (FALSE, default).
#' @param treatOther a logical stating whether elements of the \code{@@other}
#' slot should be treated as well (TRUE), or not (FALSE). See details in
#' accessor documentations (\code{\link{pop}}).
#' @param quiet a logical indicating whether warnings should be issued when
#' trying to subset components of the \code{@@other} slot (TRUE), or not (FALSE,
#' default).
#' @param \dots further arguments passed to the genlight constructor.
#' @return According to 'res.type': a list of \linkS4class{genind} object
#' (default) or a matrix of data corresponding to the 'tab' slots.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{seploc}}, \code{\link{repool}}
#' @keywords manip
#' @examples
#' 
#' \dontrun{
#' data(microbov)
#' 
#' obj <- seppop(microbov)
#' names(obj)
#' 
#' obj$Salers
#' 
#' 
#' #### example for genlight objects ####
#' x <- new("genlight", list(a=rep(1,1e3),b=rep(0,1e3),c=rep(1, 1e3)))
#' x
#' 
#' pop(x) # no population info
#' pop(x) <- c("pop1","pop2", "pop1") # set population memberships
#' pop(x)
#' seppop(x)
#' as.matrix(seppop(x)$pop1)[,1:20]
#' as.matrix(seppop(x)$pop2)[,1:20,drop=FALSE]
#' }
#' 
NULL





#' Importing data from an alignement of sequences to a genind object
#' 
#' These functions take an alignement of sequences and translate SNPs into a
#' \linkS4class{genind} object. Note that only polymorphic loci are
#' retained.\cr
#' 
#' Currently, accepted sequence formats are:\cr - DNAbin (ape package):
#' function DNAbin2genind\cr - alignment (seqinr package): function
#' alignment2genind\cr
#' 
#' 
#' @aliases DNAbin2genind alignment2genind
#' @param x an object containing aligned sequences.
#' @param pop an optional factor giving the population to which each sequence
#' belongs.
#' @param exp.char a vector of single character providing expected values; all
#' other characters will be turned to NA.
#' @param na.char a vector of single characters providing values that should be
#' considered as NA. If not NULL, this is used instead of \code{exp.char}.
#' @param polyThres the minimum frequency of a minor allele for a locus to be
#' considered as polymorphic (defaults to 0.01).
#' @return an object of the class \linkS4class{genind}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{import2genind}}, \code{\link{read.genetix}},
#' \code{\link{read.fstat}}, \code{\link{read.structure}},
#' \code{\link{read.genepop}}, \code{\link[ape]{DNAbin}},
#' \code{\link[seqinr]{as.alignment}}.
#' @keywords manip
#' @examples
#' 
#' \dontrun{
#' data(woodmouse)
#' x <- DNAbin2genind(woodmouse)
#' x
#' genind2df(x)
#' }
#' 
#' if(require(seqinr)){
#' mase.res   <- read.alignment(file=system.file("sequences/test.mase",package="seqinr"),
#' format = "mase")
#' mase.res
#' x <- alignment2genind(mase.res)
#' x
#' locNames(x) # list of polymorphic sites
#' genind2df(x)
#' 
#' ## look at Euclidean distances
#' D <- dist(truenames(x))
#' D
#' 
#' ## summarise with a PCoA
#' pco1 <- dudi.pco(D, scannf=FALSE,nf=2)
#' scatter(pco1, posi="bottomright")
#' title("Principal Coordinate Analysis\n-based on proteic distances-")
#' 
#' }
#' 
NULL





#' Web servers for adegenet
#' 
#' The function \code{adegenetServer} opens up a web page providing a simple
#' user interface for some of the functionalities implemented in adegenet.
#' These servers have been developed using the package \code{shiny}.\cr
#' 
#' Currently available servers include: \itemize{ \item \code{DAPC}: a server
#' for the Discriminant Analysis of Principal Components (see ?dapc) }
#' 
#' 
#' @aliases adegenetServer
#' @param what a character string indicating which server to start; currently
#' accepted values are: "DAPC"
#' @return The function invisibly returns NULL.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk} Caitlin Collins
#' @seealso \link{dapc}
#' @examples
#' 
#' \dontrun{
#' ## this opens a web page for DAPC
#' adegenetServer()
#' }
#' 
NULL





#' Simulated genotypes of two georeferenced populations
#' 
#' This simple data set was obtained by sampling two populations evolving in a
#' island model, simulated using Easypop (2.0.1). See \code{source} for
#' simulation details. Sample sizes were respectively 100 and 30 genotypes. The
#' genotypes were given spatial coordinates so that both populations were
#' spatially differentiated.
#' 
#' 
#' @name sim2pop
#' @docType data
#' @format \code{sim2pop} is a genind object with a matrix of xy coordinates as
#' supplementary component.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @references Balloux F (2001) Easypop (version 1.7): a computer program for
#' oppulation genetics simulations \emph{Journal of Heredity}, \bold{92}:
#' 301-302
#' @source Easypop version 2.0.1 was run with the following parameters:\cr -
#' two diploid populations, one sex, random mating\cr - 1000 individuals per
#' population\cr - proportion of migration: 0.002\cr - 20 loci\cr - mutation
#' rate: 0.0001 (KAM model)\cr - maximum of 50 allelic states\cr - 1000
#' generations (last one taken)\cr
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' data(sim2pop)
#' 
#' if(require(hierfstat)){
#' ## try and find the Fst
#' temp <- genind2hierfstat(sim2pop)
#' varcomp.glob(temp[,1],temp[,-1])
#' # Fst = 0.038
#' }
#' 
#' ## run monmonier algorithm
#' 
#' # build connection network
#' gab <- chooseCN(sim2pop@@other$xy,ask=FALSE,type=2)
#' 
#' # filter random noise 
#' pca1 <- dudi.pca(sim2pop@@tab,scale=FALSE, scannf=FALSE, nf=1)
#' 
#' # run the algorithm
#' mon1 <- monmonier(sim2pop@@other$xy,dist(pca1$l1[,1]),gab, scanthres=FALSE)
#' 
#' # graphical display
#' temp <- sim2pop@@pop
#' levels(temp) <- c(17,19)
#' temp <- as.numeric(as.character(temp))
#' plot(mon1)
#' points(sim2pop@@other$xy,pch=temp,cex=2)
#' legend("topright",leg=c("Pop A", "Pop B"),pch=c(17,19))
#' }
#' 
NULL





#' Formal class "SNPbin"
#' 
#' The class \code{SNPbin} is a formal (S4) class for storing a genotype of
#' binary SNPs in a compact way, using a bit-level coding scheme.  This storage
#' is most efficient with haploid data, where the memory taken to represent
#' data can reduced more than 50 times. However, \code{SNPbin} can be used for
#' any level of ploidy, and still remain an efficient storage mode.
#' 
#' A \code{SNPbin} object can be constructed from a vector of integers giving
#' the number of the second allele for each locus.
#' 
#' \code{SNPbin} stores a single genotype. To store multiple genotypes, use the
#' \linkS4class{genlight} class.
#' 
#' 
#' @name SNPbin-class
#' @aliases SNPbin SNPbin-class [,SNPbin-method [,SNPbin,ANY,ANY-method
#' initialize,SNPbin-method show,SNPbin-method nLoc,SNPbin-method
#' $,SNPbin-method $<-,SNPbin-method names,SNPbin-method ploidy,SNPbin-method
#' ploidy<-,SNPbin-method coerce,SNPbin,integer-method as.integer.SNPbin
#' NA.posi,SNPbin-method cbind.SNPbin c.SNPbin as,integer,SNPbin-method
#' as,numeric,SNPbin-method
#' @docType class
#' @section Objects from the class SNPbin: \code{SNPbin} objects can be created
#' by calls to \code{new("SNPbin", ...)}, where '...' can be the following
#' arguments:
#' 
#' \describe{ \item{list("snp")}{a vector of integers or numeric giving numbers
#' of copies of the second alleles for each locus. If only one unnamed argument
#' is provided to 'new', it is considered as this one.}
#' \item{list("ploidy")}{an integer indicating the ploidy of the genotype; if
#' not provided, will be guessed from the data (as the maximum from the 'snp'
#' input vector).} \item{list("label")}{an optional character string serving as
#' a label for the genotype.} }
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#' @seealso Related class:\cr - \code{\linkS4class{genlight}}, for storing
#' multiple binary SNP genotypes. \cr - \code{\linkS4class{genind}}, for
#' storing other types of genetic markers. \cr
#' @keywords classes
#' @examples
#' 
#' \dontrun{
#' #### HAPLOID EXAMPLE ####
#' ## create a genotype of 100,000 SNPs
#' dat <- sample(c(0,1,NA), 1e5, prob=c(.495, .495, .01), replace=TRUE)
#' dat[1:10]
#' x <- new("SNPbin", dat)
#' x
#' x[1:10] # subsetting
#' as.integer(x[1:10])
#' 
#' ## try a few accessors
#' ploidy(x)
#' nLoc(x)
#' head(x$snp[[1]]) # internal bit-level coding
#' 
#' ## check that conversion is OK
#' identical(as(x, "integer"),as.integer(dat)) # SHOULD BE TRUE
#' 
#' ## compare the size of the objects
#' print(object.size(dat), unit="auto")
#' print(object.size(x), unit="auto")
#' object.size(dat)/object.size(x) # EFFICIENCY OF CONVERSION
#' 
#' 
#' #### TETRAPLOID EXAMPLE ####
#' ## create a genotype of 100,000 SNPs
#' dat <- sample(c(0:4,NA), 1e5, prob=c(rep(.995/5,5), 0.005), replace=TRUE)
#' x <- new("SNPbin", dat)
#' identical(as(x, "integer"),as.integer(dat)) # MUST BE TRUE
#' 
#' ## compare the size of the objects
#' print(object.size(dat), unit="auto")
#' print(object.size(x), unit="auto")
#' object.size(dat)/object.size(x) # EFFICIENCY OF CONVERSION
#' 
#' 
#' #### c, cbind ####
#' a <- new("SNPbin", c(1,1,1,1,1))
#' b <- new("SNPbin", c(0,0,0,0,0))
#' a
#' b
#' ab <- c(a,b)
#' ab
#' identical(c(a,b),cbind(a,b))
#' as.integer(ab)
#' }
#' 
NULL





#' Analyse the position of polymorphic sites
#' 
#' These functions are used to describe the distribution of polymorphic sites
#' (SNPs) in an alignment.
#' 
#' The function \code{snpposi.plot} plots the positions and density of SNPs in
#' the alignment.
#' 
#' The function \code{snpposi.test} tests whether SNPs are randomly distributed
#' in the genome, the alternative hypothesis being that they are clustered.
#' This test is based on the distances of each SNP to the closest SNP. This
#' provides one measure of clustering for each SNP.  Different statistics can
#' be used to summarise these values (argument \code{stat}), but by default the
#' statistics used is the median.
#' 
#' \code{snpposi.plot} and \code{snpposi.test} are generic functions with
#' methods for vectors of integers or numeric (indicating SNP position), and
#' for \code{\link[ape]{DNAbin}} objects.
#' 
#' 
#' @aliases snpposi.plot snpposi.plot.integer snpposi.plot.numeric
#' snpposi.plot.DNAbin snpposi.test snpposi.test.integer snpposi.test.numeric
#' snpposi.test.DNAbin
#' @param x a vector of integers or numerics containing SNP positions, or a set
#' of aligned sequences in a \code{DNAbin} object.
#' @param genome.size an integer indicating the length of genomes.
#' @param smooth a smoothing parameter for the density estimation; smaller
#' values will give more local peaks; values have to be positive but can be
#' less than 1.
#' @param col the color to be used for the plot; ignored if codon positions are
#' represented.
#' @param alpha the alpha level to be used for transparency (density curve).
#' @param codon a logical indicating if codon position should be indicated
#' (TRUE, default) or not.
#' @param start.at an integer indicating at which base of a codon the alignment
#' starts (defaults to 1); values other than 1, 2 and 3 will be ignored.
#' @param n.sim an integer indicating the number of randomizations to be used
#' in the Monte Carlo test.
#' @param stat a function used to summarize the measure of physical proximity
#' between SNPs; by default, the median is used.
#' @param \dots further arguments to be passed to the \code{integer} method.
#' @return A Monte Carlo test of the class \code{randtest}.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}.
#' @seealso The \code{\link{fasta2DNAbin}} to read fasta alignments with
#' minimum RAM use.
#' @examples
#' 
#' if(require(ape)){
#' data(woodmouse)
#' snpposi.plot(woodmouse, codon=FALSE)
#' snpposi.plot(woodmouse)
#' 
#' \dontrun{
#' snpposi.test(c(1,3,4,5), 100)
#' snpposi.test(woodmouse)
#' }
#' }
#' 
NULL





#' Simulated data illustrating the sPCA
#' 
#' Datasets illustrating the spatial Principal Component Analysis (Jombart et
#' al. 2009).  These data were simulated using various models using Easypop
#' (2.0.1).  Spatial coordinates were defined so that different spatial
#' patterns existed in the data. The \code{spca-illus} is a list containing the
#' following \linkS4class{genind} or \linkS4class{genpop} objects:\cr - dat2A:
#' 2 patches \cr - dat2B: cline between two pop \cr - dat2C: repulsion among
#' individuals from the same gene pool \cr - dat3: cline and repulsion \cr -
#' dat4: patches and local alternance \cr
#' 
#' See "source" for a reference providing simulation details.
#' 
#' 
#' @name spcaIllus
#' @docType data
#' @format \code{spcaIllus} is list of 5 components being either genind or
#' genpop objects.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{spca}}
#' @references Jombart, T., Devillard, S., Dufour, A.-B. and Pontier, D.
#' Revealing cryptic spatial patterns in genetic variability by a new
#' multivariate method. \emph{Heredity}, \bold{101}, 92--103.
#' 
#' Balloux F (2001) Easypop (version 1.7): a computer program for oppulation
#' genetics simulations \emph{Journal of Heredity}, \bold{92}: 301-302
#' @source Jombart, T., Devillard, S., Dufour, A.-B. and Pontier, D. Revealing
#' cryptic spatial patterns in genetic variability by a new multivariate
#' method. \emph{Heredity}, \bold{101}, 92--103.
#' @keywords datasets spatial
#' @examples
#' 
#' if(require(spdep)){
#' 
#' data(spcaIllus)
#' attach(spcaIllus)
#' opar <- par(no.readonly=TRUE)
#' ## comparison PCA vs sPCA
#' 
#' # PCA
#' pca2A <- dudi.pca(dat2A$tab,center=TRUE,scale=FALSE,scannf=FALSE)
#' pca2B <- dudi.pca(dat2B$tab,center=TRUE,scale=FALSE,scannf=FALSE)
#' pca2C <- dudi.pca(dat2C$tab,center=TRUE,scale=FALSE,scannf=FALSE)
#' pca3 <- dudi.pca(dat3$tab,center=TRUE,scale=FALSE,scannf=FALSE,nf=2)
#' pca4 <- dudi.pca(dat4$tab,center=TRUE,scale=FALSE,scannf=FALSE,nf=2)
#' 
#' # sPCA
#' spca2A <-spca(dat2A,xy=dat2A$other$xy,ask=FALSE,type=1,
#' plot=FALSE,scannf=FALSE,nfposi=1,nfnega=0)
#' 
#' spca2B <- spca(dat2B,xy=dat2B$other$xy,ask=FALSE,type=1,
#' plot=FALSE,scannf=FALSE,nfposi=1,nfnega=0)
#' 
#' spca2C <- spca(dat2C,xy=dat2C$other$xy,ask=FALSE,
#' type=1,plot=FALSE,scannf=FALSE,nfposi=0,nfnega=1)
#' 
#' spca3 <- spca(dat3,xy=dat3$other$xy,ask=FALSE,
#' type=1,plot=FALSE,scannf=FALSE,nfposi=1,nfnega=1)
#' 
#' spca4 <- spca(dat4,xy=dat4$other$xy,ask=FALSE,
#' type=1,plot=FALSE,scannf=FALSE,nfposi=1,nfnega=1)
#' 
#' # an auxiliary function for graphics
#' plotaux <- function(x,analysis,axis=1,lab=NULL,...){
#' neig <- NULL
#' if(inherits(analysis,"spca")) neig <- nb2neig(analysis$lw$neighbours)
#' xrange <- range(x$other$xy[,1])
#' xlim <- xrange + c(-diff(xrange)*.1 , diff(xrange)*.45)
#' yrange <- range(x$other$xy[,2])
#' ylim <- yrange + c(-diff(yrange)*.45 , diff(yrange)*.1)
#' 
#' s.value(x$other$xy,analysis$li[,axis],include.ori=FALSE,addaxes=FALSE,
#' cgrid=0,grid=FALSE,neig=neig,cleg=0,xlim=xlim,ylim=ylim,...)
#' 
#' par(mar=rep(.1,4))
#' if(is.null(lab)) lab = gsub("[P]","",x$pop)
#' text(x$other$xy, lab=lab, col="blue", cex=1.2, font=2)
#' add.scatter({barplot(analysis$eig,col="grey");box();
#' title("Eigenvalues",line=-1)},posi="bottomright",ratio=.3)
#' }
#' 
#' # plots
#' plotaux(dat2A,pca2A,sub="dat2A - PCA",pos="bottomleft",csub=2)
#' plotaux(dat2A,spca2A,sub="dat2A - sPCA glob1",pos="bottomleft",csub=2)
#' 
#' plotaux(dat2B,pca2B,sub="dat2B - PCA",pos="bottomleft",csub=2)
#' plotaux(dat2B,spca2B,sub="dat2B - sPCA glob1",pos="bottomleft",csub=2)
#' 
#' plotaux(dat2C,pca2C,sub="dat2C - PCA",pos="bottomleft",csub=2)
#' plotaux(dat2C,spca2C,sub="dat2C - sPCA loc1",pos="bottomleft",csub=2,axis=2)
#' 
#' par(mfrow=c(2,2))
#' plotaux(dat3,pca3,sub="dat3 - PCA axis1",pos="bottomleft",csub=2)
#' plotaux(dat3,spca3,sub="dat3 - sPCA glob1",pos="bottomleft",csub=2)
#' plotaux(dat3,pca3,sub="dat3 - PCA axis2",pos="bottomleft",csub=2,axis=2)
#' plotaux(dat3,spca3,sub="dat3 - sPCA loc1",pos="bottomleft",csub=2,axis=2)
#' 
#' plotaux(dat4,pca4,lab=dat4$other$sup.pop,sub="dat4 - PCA axis1",
#' pos="bottomleft",csub=2)
#' plotaux(dat4,spca4,lab=dat4$other$sup.pop,sub="dat4 - sPCA glob1",
#' pos="bottomleft",csub=2)
#' plotaux(dat4,pca4,lab=dat4$other$sup.pop,sub="dat4 - PCA axis2",
#' pos="bottomleft",csub=2,axis=2)
#' plotaux(dat4,spca4,lab=dat4$other$sup.pop,sub="dat4 - sPCA loc1",
#' pos="bottomleft",csub=2,axis=2)
#' 
#' # color plot
#' par(opar)
#' colorplot(spca3, cex=4, main="colorplot sPCA dat3")
#' text(spca3$xy[,1], spca3$xy[,2], dat3$pop)
#' 
#' colorplot(spca4, cex=4, main="colorplot sPCA dat4")
#' text(spca4$xy[,1], spca4$xy[,2], dat4$other$sup.pop)
#' 
#' # detach data
#' detach(spcaIllus)
#' }
#' 
NULL





#' Virtual classes for adegenet
#' 
#' These virtual classes are only for internal use in adegenet
#' 
#' 
#' @name virtualClasses
#' @aliases indInfo-class popInfo-class gen-class callOrNULL-class
#' charOrNULL-class factorOrNULL-class intOrNum-class listOrNULL-class
#' intOrNULL-class
#' @docType class
#' @section Objects from the Class: A virtual Class: No objects may be created
#' from it.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords classes
NULL





#' Website and tutorials for adegenet
#' 
#' These functions allow to access documentation for adegenet available
#' online.\cr
#' 
#' These functions include: \itemize{ \item \code{adegenetWeb}: opens the
#' adegenet website in a web navigator.  \item \code{adegenetTutorial}: opens a
#' specified tutorial.  }
#' 
#' 
#' @name Website and tutorials
#' @aliases adegenetWeb adegenetTutorial adegenetTutorials
#' @docType methods
#' @param which a character string indicating the type of tutorial to open
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @examples
#' 
#' \dontrun{
#' ## this opens the adegenet website
#' adegenetWeb()
#' adegenetTutorial("dapc")
#' }
#' 
NULL





#' Cross-validation for Discriminant Analysis of Principal Components (DAPC)
#' 
#' The function \code{xvalDapc} performs stratified cross-validation of DAPC
#' using varying numbers of PCs (and keeping the number of discriminant
#' functions fixed); \code{xvalDapc} is a generic with methods for
#' \code{data.frame} and \code{matrix}.\cr
#' 
#' The Discriminant Analysis of Principal Components (DAPC) relies on dimension
#' reduction of the data using PCA followed by a linear discriminant analysis.
#' How many PCA axes to retain is often a non-trivial question. Cross
#' validation provides an objective way to decide how many axes to retain:
#' different numbers are tried and the quality of the corresponding DAPC is
#' assessed by cross-validation: DAPC is performed on a training set, typically
#' made of 90\% of the observations (comprising 90\% of the observations in
#' each subpopulation) , and then used to predict the groups of the 10\% of
#' remaining observations.  The current method uses the average prediction
#' success per group (result="groupMean"), or the overall prediction success
#' (result="overall"). The number of PCs associated with the lowest Mean
#' Squared Error is then retained in the DAPC.
#' 
#' @aliases xvalDapc xvalDapc.data.frame xvalDapc.matrix
#' @param x \code{a data.frame} or a \code{matrix} used as input of DAPC.
#' @param grp a \code{factor} indicating the group membership of individuals.
#' @param n.pca.max maximum number of PCA components to retain.
#' @param n.da an \code{integer} indicating the number of axes retained in the
#' Discriminant Analysis step. If \code{NULL}, n.da defaults to 1 less than the
#' number of groups.
#' @param training.set the proportion of data (individuals) to be used for the
#' training set; defaults to 0.9 if all groups have >= 10 members; otherwise,
#' training.set scales automatically to the largest proportion that still
#' ensures all groups will be present in both training and validation sets.
#' @param result a character string; "groupMean" for group-wise assignment
#' sucess, or "overall" for an overall mean assignment success; see details.
#' @param center a \code{logical} indicating whether variables should be
#' centred to mean 0 (TRUE, default) or not (FALSE). Always TRUE for
#' \linkS4class{genind} objects.
#' @param scale a \code{logical} indicating whether variables should be scaled
#' (TRUE) or not (FALSE, default). Scaling consists in dividing variables by
#' their (estimated) standard deviation to account for trivial differences in
#' variances.
#' @param n.pca an \code{integer} vector indicating the number of different
#' number of PCA axes to be retained for the cross validation; if \code{NULL},
#' this will be dertermined automatically.
#' @param n.rep the number of replicates to be carried out at each level of PC
#' retention; defaults to 30.
#' @param xval.plot a logical indicating whether a plot of the cross-validation
#' results should be generated.
#' @param \dots further arguments to be passed to other methods.
#' @return A \code{list} containing seven items, and a \code{plot} of the
#' results.  The first is a \code{data.frame} with two columns, the first
#' giving the number of PCs of PCA retained in the corresponding DAPC, and the
#' second giving the proportion of successful group assignment for each
#' replicate.  The second item gives the mean and confidence interval for
#' random chance.  The third gives the mean successful assignment at each level
#' of PC retention.  The fourth indicates which number of PCs is associated
#' with the highest mean success.  The fifth gives the Root Mean Squared Error
#' at each level of PC retention.  The sixth indicates which number of PCs is
#' associated with the lowest MSE.  The seventh item contains the DAPC carried
#' out with the optimal number of PCs, determined with reference to MSE.
#' 
#' If \code{xval.plot=TRUE} a scatterplot of the results of cross-validation
#' will be displayed.
#' @author Caitlin Collins \email{caitlin.collins12@@imperial.ac.uk}, Thibaut
#' Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{dapc}}
#' @references Jombart T, Devillard S and Balloux F (2010) Discriminant
#' analysis of principal components: a new method for the analysis of
#' genetically structured populations. BMC Genetics11:94.
#' doi:10.1186/1471-2156-11-94
#' @keywords multivariate
#' @examples
#' 
#' \dontrun{
#' ## CROSS-VALIDATION ##
#' data(sim2pop)
#' xval <- xvalDapc(sim2pop@@tab, pop(sim2pop), n.pca.max=100, n.rep=3)
#' xval
#' 
#' }
#' 
NULL



