## ###############
## ## simOutbreak
## ###############
## simOutbreak <- function(R0, infec.curve, n.hosts=200, duration=50,
##                         seq.length=1e4, mu.transi=1e-4, mu.transv=mu.transi/2,
##                         tree=TRUE){

##     ## CHECKS ##
##     if(!require(ape)) stop("The ape package is required.")


##     ## HANDLE ARGUMENTS ##
##     ## normalize gen.time
##     infec.curve <- infec.curve/sum(infec.curve)
##     infec.curve <- c(infec.curve, rep(0, duration)) # make sure dates go all the way
##     t.clear <- which(diff(infec.curve<1e-10)==1) # time at which infection is cleared

##     ## GENETIC FUNCTIONS ##
##     NUCL <- as.DNAbin(c("a","t","c","g"))
##     TRANSISET <- list('a'=as.DNAbin('g'), 'g'=as.DNAbin('a'), 'c'=as.DNAbin('t'), 't'=as.DNAbin('c'))
##     TRANSVSET <- list('a'=as.DNAbin(c('c','t')), 'g'=as.DNAbin(c('c','t')), 'c'=as.DNAbin(c('a','g')), 't'=as.DNAbin(c('a','g')))


##     ## AUXILIARY FUNCTIONS ##
##     ## generate sequence from scratch
##     seq.gen <- function(){
##         ##res <- list(sample(NUCL, size=seq.length, replace=TRUE)) # DNAbin are no longer lists by default
##         res <- sample(NUCL, size=seq.length, replace=TRUE)
##         class(res) <- "DNAbin"
##         return(res)
##     }

##     ## create substitutions for defined SNPs - no longer used
##     substi <- function(snp){
##         res <- sapply(1:length(snp), function(i) sample(setdiff(NUCL,snp[i]),1)) # ! sapply does not work on DNAbin vectors directly
##         class(res) <- "DNAbin"
##         return(res)
##     }

##     ## create transitions for defined SNPs
##     transi <- function(snp){
##         res <- unlist(TRANSISET[as.character(snp)])
##         class(res) <- "DNAbin"
##         return(res)
##     }

##     ## create transversions for defined SNPs
##     transv <- function(snp){
##         res <- sapply(TRANSVSET[as.character(snp)],sample,1)
##         class(res) <- "DNAbin"
##         return(res)
##     }

##     ## duplicate a sequence (including possible mutations)
##     seq.dupli <- function(seq, T){ # T is the number of time units between ancestor and decendent
##         ## transitions ##
##         n.transi <- rbinom(n=1, size=seq.length*T, prob=mu.transi) # total number of transitions
##         if(n.transi>0) {
##             idx <- sample(1:seq.length, size=n.transi, replace=FALSE)
##             seq[idx] <- transi(seq[idx])
##         }

##         ## transversions ##
##         n.transv <- rbinom(n=1, size=seq.length*T, prob=mu.transv) # total number of transitions
##         if(n.transv>0) {
##             idx <- sample(1:seq.length, size=n.transv, replace=FALSE)
##             seq[idx] <- transv(seq[idx])
##         }
##         return(seq)
##     }


##     ## MAIN FUNCTION ##
##     ## initialize results ##
##     dynam <- data.frame(nsus=integer(duration+1), ninf=integer(duration+1), nrec=integer(duration+1))
##     rownames(dynam) <- 0:duration
##     res <- list(n=1, dna=NULL, dates=NULL, id=NULL, ances=NULL, dynam=dynam)
##     res$dynam$nsus[1] <- n.hosts-1
##     res$dynam$ninf[1] <- 1
##     res$dates[1] <- 0
##     res$ances <- NA
##     res$dna <- matrix(seq.gen(),nrow=1)
##     class(res$dna) <- "DNAbin"

##     ## run outbreak ##
##     for(t in 1:duration){
##         ## individual force of infection
##         indivForce <- infec.curve[t-res$dates+1]

##         ## global force of infection (R0 \sum_j I_t^j / N)
##         globForce <- sum(indivForce)*R0/n.hosts

##         ## number of new infections
##         nbNewInf <- rbinom(1, size=res$dynam$nsus[t], prob=globForce)

##         ## dates of new infections
##         if(nbNewInf>0){
##             res$dates <- c(res$dates, rep(t,nbNewInf))

##             ## ancestries of the new infections
##             temp <- as.vector(rmultinom(1, size=nbNewInf, prob=indivForce))
##             newAnces <- rep(which(temp>0), temp[which(temp>0)])
##             res$ances <- c(res$ances,newAnces)

##             ## dna sequences of the new infections
##             newSeq <- t(sapply(newAnces, function(i) seq.dupli(res$dna[i,], t-res$dates[i])))
##             res$dna <- rbind(res$dna, newSeq)
##         }

##         ## update nb of infected, recovered, etc.
##         res$dynam$nrec[t+1] <- sum(res$dates>=t.clear)
##         res$dynam$ninf[t+1] <- sum(res$dates>=0 & res$dates < t.clear)
##         res$dynam$nsus[t+1] <- res$dynam$nsus[t] - nbNewInf
##     } # end for


##     ## SHAPE AND RETURN OUTPUT ##
##     res$n <- nrow(res$dna)
##     res$id <- 1:res$n
##     res$nmut <- sapply(1:res$n, function(i) dist.dna(res$dna[c(res$id[i],res$ances[i]),], model="raw"))*ncol(res$dna)
##     res$call <- match.call()
##     if(tree){
##         res$tree <- fastme.ols(dist.dna(res$dna, model="TN93"))
##         res$tree <- root(res$tree,"1")
##     }
##     class(res) <- "simOutbreak"
##     return(res)

## } # end simOutbreak








## ##################
## ## print.simOutbreak
## ##################
## print.simOutbreak <- function(x, ...){

##     cat("\t\n=========================")
##     cat("\t\n=   simulated outbreak  =")
##     cat("\t\n=  (simOutbreak object) =")
##     cat("\t\n=========================\n")

##     cat("\nSize :", x$n,"cases (out of", x$dynam$nsus[1],"susceptible hosts)")
##     cat("\nGenome length :", ncol(x$dna),"nucleotids")
##     cat("\nDate range :", min(x$dates),"-",max(x$dates))

##     cat("\nContent:\n")
##     print(names(x))

##     return(NULL)
## } # end print.simOutbreak







## ##############
## ## [.simOutbreak
## ##############
## "[.simOutbreak" <- function(x,i,j,drop=FALSE){
##     res <- x
##     res$dna <- res$dna[i,,drop=FALSE]
##     res$id <- res$id[i]
##     res$ances <- res$ances[i]
##     res$ances[!res$ances %in% res$id] <- NA
##     res$dates <- res$dates[i]
##     res$n <- nrow(res$dna)

##     return(res)
## }





## ##################
## ## labels.simOutbreak
## ##################
## labels.simOutbreak <- function(object, ...){
##     return(object$id)
## }






## #########################
## ## as.igraph.simOutbreak
## #########################
## as.igraph.simOutbreak <- function(x, ...){
##     if(!require(igraph)) stop("package igraph is required for this operation")
##     if(!require(ape)) stop("package ape is required for this operation")

##     ## GET DAG ##
##     from <- x$ances
##     to <- x$id
##     isNotNA <- !is.na(from) & !is.na(to)
##     dat <- data.frame(from,to,stringsAsFactors=FALSE)[isNotNA,,drop=FALSE]
##     vnames <- as.character(unique(unlist(dat)))
##     out <- graph.data.frame(dat, directed=TRUE, vertices=data.frame(names=vnames, dates=x$dates[vnames]))

##     ## SET WEIGHTS ##
##     D <- as.matrix(dist.dna(x$dna,model="raw")*ncol(x$dna))
##     temp <- mapply(function(i,j) return(D[i,j]), as.integer(from), as.integer(to))
##     E(out)$weight <- temp[isNotNA]

##     ## SET ARROW WIDTH ##
##     temp <- max(E(out)$weight) - E(out)$weight
##     temp <- temp/max(temp) * 4
##     E(out)$width <- round(temp)+1

##     return(out)
## }







## ####################
## ## plot.simOutbreak
## ####################
## plot.simOutbreak <- function(x, y=NULL, cex=1, col=num2col(x$dates), label=x$id,
##                              edge.col=num2col(x$nmut[-1], col.pal=seasun), lwd=1, ...){
##     if(!require(igraph)) stop("package igraph is required for this operation")
##     if(!require(ape)) stop("package ape is required for this operation")
##     plot(as.igraph(x), vertex.size=15*cex, vertex.color=col, vertex.label=label,
##          vertex.label.cex=cex, edge.color=edge.col, edge.width=lwd, ...)
## } # end plot.simOutbreak




















## ## #####################
## ## ## seqTrack.simOutbreak
## ## #####################
## ## seqTrack.simOutbreak <- function(x, best=c("min","max"), prox.mat=NULL, ...){
## ##     myX <- dist.dna(x$dna, model="raw")
## ##     x.names <- labels(x)
## ##     x.dates <- as.POSIXct(x)
## ##     seq.length <- ncol(x$dna)
## ##     myX <- myX * seq.length
## ##     myX <- as.matrix(myX)
## ##     prevCall <- as.list(x$call)
## ##     if(is.null(prevCall$mu)){
## ##         mu0 <- 0.0001
## ##     } else {
## ##         mu0 <- eval(prevCall$mu)
## ##     }
## ##     res <- seqTrack(myX, x.names=x.names, x.dates=x.dates, best=best, prox.mat=prox.mat,...)
## ##     return(res)
## ## }






## ## ########################
## ## ## as.seqTrack.simOutbreak
## ## ########################
## ## as.seqTrack.simOutbreak <- function(x){
## ##     ## x.ori <- x
## ##     ## x <- na.omit(x)
## ##     toSetToNA <- x$dates==min(x$dates)
## ##     res <- list()
## ##     res$id <- labels(x)
## ##     res <- as.data.frame(res)
## ##     res$ances <- x$ances
## ##     res$ances[toSetToNA] <- NA
## ##     res$weight <- 1 # ??? have to recompute that...
## ##     res$weight[toSetToNA] <- NA
## ##     res$date <- as.POSIXct(x)[labels(x)]
## ##     res$ances.date <- as.POSIXct(x)[x$ances]

## ##     ## set results as indices rather than labels
## ##     res$ances <- match(res$ances, res$id)
## ##     res$id <- 1:length(res$id)

## ##     ## SET CLASS
## ##     class(res) <- c("seqTrack", "data.frame")

## ##     return(res)
## ## }








## ## ###################
## ## ## sample.simOutbreak
## ## ###################
## ## sample.simOutbreak <- function(x, n){
## ## ##sample.simOutbreak <- function(x, n, rDate=.rTimeSeq, arg.rDate=NULL){
## ##     ## EXTRACT THE SAMPLE ##
## ##     res <- x[sample(1:nrow(x$dna), n, replace=FALSE)]


## ##     ## RETRIEVE SOME PARAMETERS FROM HAPLOSIM CALL
## ##     prevCall <- as.list(x$call)
## ##     if(!is.null(prevCall$mu)){
## ##         mu0 <- eval(prevCall$mu)
## ##     } else {
## ##         mu0 <- 1e-04
## ##     }

## ##     if(!is.null(prevCall$dna.length)){
## ##         L <- eval(prevCall$dna.length)
## ##     } else {
## ##         L <- 1000
## ##     }

## ##     ## truedates <- res$dates
## ##     ## daterange <- diff(range(res$dates,na.rm=TRUE))

## ##     ## if(identical(rDate,.rTimeSeq)){
## ##     ##     sampdates <- .rTimeSeq(n=length(truedates), mu=mu0, L=L, maxNbDays=daterange/2)
## ##     ## } else{
## ##     ##     arg.rDate$n <- n
## ##     ##     sampdates <- do.call(rDate, arg.rDate)
## ##     ## }
## ##     ## sampdates <- truedates + abs(sampdates)

## ##     return(res)
## ## } # end sample.simOutbreak














