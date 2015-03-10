\encoding{UTF-8}
\name{microbov}
\alias{microbov}
\docType{data}
\title{Microsatellites genotypes of 15 cattle breeds}
\description{
This data set gives the genotypes of 704 cattle individuals for 30
microsatellites recommended by the FAO. The individuals are divided into
two countries (Afric, France), two species (Bos taurus, Bos indicus) and
15 breeds. Individuals were chosen in order to avoid pseudoreplication
according to their exact genealogy.
}
\usage{data(microbov)}
\format{
    \code{microbov} is a genind object with 3 supplementary components:
    \describe{
        \item{coun}{a factor giving the country of each individual (AF:
	  Afric; FR: France).}
        \item{breed}{a factor giving the breed of each individual.}
        \item{spe}{is a factor giving the species of each individual
	  (BT: Bos taurus; BI: Bos indicus).}
    }
}
\source{
Data prepared by Katayoun Moazami-Goudarzi and Denis Lalo\"e (INRA,
Jouy-en-Josas, France)
}
\references{
  Lalo\"e D., Jombart T., Dufour A.-B. and Moazami-Goudarzi K. (2007)
  Consensus genetic structuring and typological value of markers using
  Multiple Co-Inertia Analysis. \emph{Genetics Selection Evolution}.
  \bold{39}: 545--567.
}
\examples{
\dontrun{
data(microbov)
microbov
summary(microbov)

# make Y, a genpop object
Y <- genind2genpop(microbov)

# make allelic frequency table
temp <- makefreq(Y,missing="mean")
X <- temp$tab
nsamp <- temp$nobs

# perform 1 PCA per marker 

kX <- ktab.data.frame(data.frame(X),Y@loc.nall)

kpca <- list()
for(i in 1:30) {kpca[[i]] <- dudi.pca(kX[[i]],scannf=FALSE,nf=2,center=TRUE,scale=FALSE)}


sel <- sample(1:30,4)
col = rep('red',15)
col[c(2,10)] = 'darkred'
col[c(4,12,14)] = 'deepskyblue4'
col[c(8,15)] = 'darkblue'

# display %PCA
par(mfrow=c(2,2))
for(i in sel) {
s.multinom(kpca[[i]]$c1,kX[[i]],n.sample=nsamp[,i],coulrow=col,sub=Y@loc.names[i])
add.scatter.eig(kpca[[i]]$eig,3,xax=1,yax=2,posi="top")
}

# perform a Multiple Coinertia Analysis
kXcent <- kX
for(i in 1:30) kXcent[[i]] <- as.data.frame(scalewt(kX[[i]],center=TRUE,scale=FALSE))
mcoa1 <- mcoa(kXcent,scannf=FALSE,nf=3, option="uniform")

# coordinated %PCA
mcoa.axes <- split(mcoa1$axis,Y@loc.fac)
mcoa.coord <- split(mcoa1$Tli,mcoa1$TL[,1])
var.coord <- lapply(mcoa.coord,function(e) apply(e,2,var))

par(mfrow=c(2,2))
for(i in sel) {
s.multinom(mcoa.axes[[i]][,1:2],kX[[i]],n.sample=nsamp[,i],coulrow=col,sub=Y@loc.names[i])
add.scatter.eig(var.coord[[i]],2,xax=1,yax=2,posi="top")
}

# reference typology
par(mfrow=c(1,1))
s.label(mcoa1$SynVar,lab=microbov@pop.names,sub="Reference typology",csub=1.5)
add.scatter.eig(mcoa1$pseudoeig,nf=3,xax=1,yax=2,posi="top")

# typologial values
tv <- mcoa1$cov2
tv <- apply(tv,2,function(c) c/sum(c))*100
rownames(tv) <- Y@loc.names
tv <- tv[order(Y@loc.names),]

par(mfrow=c(3,1),mar=c(5,3,3,4),las=3)
for(i in 1:3){
barplot(round(tv[,i],3),ylim=c(0,12),yaxt="n",main=paste("Typological value -
structure",i))
axis(side=2,at=seq(0,12,by=2),labels=paste(seq(0,12,by=2),"\%"),cex=3)
abline(h=seq(0,12,by=2),col="grey",lty=2)
}
}
}
\keyword{datasets}
