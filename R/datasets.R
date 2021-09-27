




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
#' kX <- ktab.data.frame(data.frame(X),Y@@loc.n.all)
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
#' s.multinom(kpca[[i]]$c1,kX[[i]],n.sample=nsamp[,i],coulrow=col,sub=locNames(Y)[i])
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
#' s.multinom(mcoa.axes[[i]][,1:2],kX[[i]],n.sample=nsamp[,i],coulrow=col,sub=locNames(Y)[i])
#' add.scatter.eig(var.coord[[i]],2,xax=1,yax=2,posi="top")
#' }
#'
#' # reference typology
#' par(mfrow=c(1,1))
#' s.label(mcoa1$SynVar,lab=popNames(microbov),sub="Reference typology",csub=1.5)
#' add.scatter.eig(mcoa1$pseudoeig,nf=3,xax=1,yax=2,posi="top")
#'
#' # typologial values
#' tv <- mcoa1$cov2
#' tv <- apply(tv,2,function(c) c/sum(c))*100
#' rownames(tv) <- locNames(Y)
#' tv <- tv[order(locNames(Y)),]
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
#' colonies as a supplementary components (@@xy).
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
#' barplot(x$loc.n.all,ylab="Alleles numbers",main="Alleles numbers
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
#' @encoding utf-8
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
#'
#' data(rupica)
#' rupica
#'
#'
#' \dontrun{
#' required_packages <- require(adehabitat) &&
#'   require(adespatial) &&
#'   require(spdep)
#' if (required_packages) {
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
#' required_packages <- require(adespatial) && require(spdep)
#' if (required_packages) {
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



#' Toy hybrid dataset
#' @name hybridtoy
#' @aliases hybridtoy
#' @docType data
#' @format a \linkS4class{genind} object
#' @author Data simulated by  Marie-Pauline Beugin. Example by Thibaut Jombart.
#'
#' @examples
#' data(hybridtoy)
#' x <- hybridtoy
#' pca1 <- dudi.pca(tab(x), scannf=FALSE, scale=FALSE)
#' s.class(pca1$li, pop(x))
#'
#' if(require(ggplot2)) {
#' p <- ggplot(pca1$li, aes(x=Axis1)) +
#'     geom_density(aes(fill=pop(x)), alpha=.4, adjust=1) +
#'     geom_point(aes(y=0, color=pop(x)), pch="|", size=10, alpha=.5)
#' p
#' }
#'
#' ## kmeans
#' km <- find.clusters(x, n.pca=10, n.clust=2)
#' table(pop(x), km$grp)
#'
#' ## dapc
#' dapc1 <- dapc(x, pop=km$grp, n.pca=10, n.da=1)
#' scatter(dapc1)
#' scatter(dapc1, grp=pop(x))
#' compoplot(dapc1, col.pal=spectral, n.col=2)
#'
#' ## ML-EM with hybrids
#' res <- snapclust(x, k=2, hybrids=TRUE, detailed=TRUE)
#' compoplot(res, n.col=3)
#' table(res$group, pop(x))
#'
NULL




#' Microsatellites genotypes of 781 swallowtail butterflies from 40 populations in
#' Alberta and British Columbia, Canada
#'
#' This data set gives the genotypes of 781 swallowtail butterflies 
#' (\emph{Papilio machaon} species group) for 10 microsatellites markers.
#' The individuals are divided into 40 populations.
#'
#'
#' @name swallowtails
#' @docType data
#' @format \code{swallowtails} is a genind object containing 781 individuals, 
#' 10 microsatellite markers, and 40 populations.
#' @references Dupuis, J.R. & Sperling, F.A.H. Hybrid dynamics in a species
#' group of swallowtail butterflies. \emph{Journal of Evolutionary Biology}, 
#' \bold{10}, 1932--1951.
#' @source Julian Dupuis (University of Hawaii, USA)
#' @keywords datasets
#' @examples
#'
#' \dontrun{
#' data(swallowtails)
#' swallowtails
#' 
#' # conducting a DAPC (n.pca determined using xvalDapc, see ??xvalDapc)
#' 
#' dapc1 <- dapc(swallowtails, n.pca=40, n.da=200)
#' 
#' # read in swallowtails_loc.csv, which contains "key", "lat", and "lon"
#' # columns with column headers (this example contains additional columns
#' # containing species identifications, locality descriptions, and COI
#' # haplotype clades)
#' 
#' input_locs <- system.file("files/swallowtails_loc.csv", package = "adegenet")
#' loc <- read.csv(input_locs, header = TRUE)
#' 
#' # generate mvmapper input file, automatically write the output to a csv, and
#' # name the output csv "mvMapper_Data.csv"
#' 
#' out <- export_to_mvmapper(dapc1, loc, write_file = TRUE, out_file = "mvMapper_Data.csv")
#' }
NULL

