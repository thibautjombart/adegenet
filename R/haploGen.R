############
## haploGen
############
##
## N: number of sequences to simulate
## mu: mutation rate per nucleotid per generation
## Tmax: periode of time to simulate
## mean.gen.time, sd.gen.time: average time for transmission and its standard deviation (normal dist)
## mean.repro, sd.repro: average number of transmissions and its standard deviation (normal dist)
##
haploGen <- function(seq.length=1e4, mu.transi=1e-4, mu.transv=mu.transi/2, t.max=20,
                     gen.time=function(){1+rpois(1,0.5)},
                     repro=function(){rpois(1,1.5)}, max.nb.haplo=200,
                     geo.sim=FALSE, grid.size=10, lambda.xy=0.5,
                     mat.connect=NULL,
                     ini.n=1, ini.xy=NULL){

    ## CHECKS ##
    ## if(!require(ape)) stop("The ape package is required.")


    ## HANDLE ARGUMENTS ##
    ## if numeric value, make it a function
    if(is.numeric(gen.time)){
        gen.time.val <- gen.time[1]
        gen.time <- function(){return(gen.time.val)}
    }

    ## if numeric value, make it a function
    if(is.numeric(repro)){
        repro.val <- repro[1]
        repro <- function(){return(repro.val)}
    }



    ## GENERAL VARIABLES ##
    NUCL <- as.DNAbin(c("a","t","c","g"))
    TRANSISET <- list('a'=as.DNAbin('g'), 'g'=as.DNAbin('a'), 'c'=as.DNAbin('t'), 't'=as.DNAbin('c'))
    TRANSVSET <- list('a'=as.DNAbin(c('c','t')), 'g'=as.DNAbin(c('c','t')), 'c'=as.DNAbin(c('a','g')), 't'=as.DNAbin(c('a','g')))
    res <- list(seq=as.matrix(as.DNAbin(character(0))), dates=integer(), ances=character())
    toExpand <- logical()
    myGrid <- matrix(1:grid.size^2, ncol=grid.size, nrow=grid.size)


    ## AUXILIARY FUNCTIONS ##
    ## generate sequence from scratch
    seq.gen <- function(){
        ##res <- list(sample(NUCL, size=seq.length, replace=TRUE)) # DNAbin are no longer lists by default
        res <- sample(NUCL, size=seq.length, replace=TRUE)
        class(res) <- "DNAbin"
        return(res)
    }

    ## create substitutions for defined SNPs - no longer used
    substi <- function(snp){
        res <- sapply(1:length(snp), function(i) sample(setdiff(NUCL,snp[i]),1)) # ! sapply does not work on DNAbin vectors directly
        class(res) <- "DNAbin"
        return(res)
    }

    ## create transitions for defined SNPs
    transi <- function(snp){
        res <- unlist(TRANSISET[as.character(snp)])
        class(res) <- "DNAbin"
        return(res)
    }

    ## create transversions for defined SNPs
    transv <- function(snp){
        res <- sapply(TRANSVSET[as.character(snp)],sample,1)
        class(res) <- "DNAbin"
        return(res)
    }

    ## duplicate a sequence (including possible mutations)
    seq.dupli <- function(seq, T){ # T is the number of time units between ancestor and decendent
        ## transitions ##
        n.transi <- rbinom(n=1, size=seq.length*T, prob=mu.transi) # total number of transitions
        if(n.transi>0) {
            idx <- sample(1:seq.length, size=n.transi, replace=FALSE)
            seq[idx] <- transi(seq[idx])
        }

        ## transversions ##
        n.transv <- rbinom(n=1, size=seq.length*T, prob=mu.transv) # total number of transitions
        if(n.transv>0) {
            idx <- sample(1:seq.length, size=n.transv, replace=FALSE)
            seq[idx] <- transv(seq[idx])
        }
        return(seq)
    }

    ## what is the name of the new sequences?
    seqname.gen <- function(nb.new.seq){
        res <- max(as.integer(rownames(res$seq)), 0) + 1:nb.new.seq
        return(as.character(res))
    }

    ## how many days before duplication occurs ?
    time.dupli <- function(){
        ##res <- round(rnorm(1, mean=mean.gen.time, sd=sd.gen.time))
        res <- round(gen.time()) # force integers
        res[res<0] <- 0
        return(res)
    }

    ## when duplication occurs?
    date.dupli <- function(curDate){
        res <- curDate + time.dupli()
        return(res)
    }

    ## how many duplication/transmission occur?
    nb.desc <- function(){
        ##res <- round(rnorm(1, mean=mean.repro, sd=sd.repro))
        res <- repro()
        res[res<0] <- 0
        return(res)
    }

    ## where does an haplotype emerges in the first place?
    xy.gen <- function(){
        return(sample(1:grid.size, size=2, replace=TRUE))
    }

    ## where does a transmission occur (destination)?
    if(is.null(mat.connect)){ # use universal lambda param
        xy.dupli <- function(cur.xy, nbLoc){
            mvt <- rpois(2*nbLoc, lambda.xy) * sample(c(-1,1), size=2*nbLoc, replace=TRUE)
            res <- t(matrix(mvt, nrow=2) + as.vector(cur.xy))
            res[res < 1] <- 1
            res[res > grid.size] <- grid.size
            return(res)
        }
    } else { # use location-dependent proba of dispersal between locations
        if(any(mat.connect < 0)) stop("Negative values in mat.connect (probabilities expected!)")
        mat.connect <- prop.table(mat.connect,1)
        xy.dupli <- function(cur.xy, nbLoc){
            idxAncesLoc <- myGrid[cur.xy[1], cur.xy[2]]
            newLoc <- sample(1:grid.size^2, size=nbLoc, prob=mat.connect[idxAncesLoc,], replace=TRUE) # get new locations
            res <- cbind(row(myGrid)[newLoc] , col(myGrid)[newLoc]) # get coords of new locations
            return(res)
        }
    }


    ## check result size and resize it if needed
    resize.result <- function(){
        curSize <- length(res$dates)
        if(curSize <= max.nb.haplo) return(NULL)
        toKeep <- rep(FALSE, curSize)
        toKeep[sample(1:curSize, size=max.nb.haplo, replace=FALSE)] <- TRUE
        removed.strains <- rownames(res$seq)[!toKeep]
        res$seq <<- res$seq[toKeep,]
        res$dates <<- res$dates[toKeep]
        res$ances <<- res$ances[toKeep]
        toExpand <<- toExpand[toKeep]
        temp <- as.character(res$ances) %in% removed.strains
        if(any(temp)) {
            res$ances[temp] <<- NA
        }

        return(NULL)
    }

    ## check result size and resize it if needed - spatial version
    resize.result.xy <- function(){
        curSize <- length(res$dates)
        if(curSize <= max.nb.haplo) return(NULL)
        toKeep <- rep(FALSE, curSize)
        toKeep[sample(1:curSize, size=max.nb.haplo, replace=FALSE)] <- TRUE
        removed.strains <- rownames(res$seq)[!toKeep]
        res$seq <<- res$seq[toKeep,]
        res$dates <<- res$dates[toKeep]
        res$ances <<- res$ances[toKeep]
        res$xy <<- res$xy[toKeep,,drop=FALSE]
        toExpand <<- toExpand[toKeep]
        temp <- as.character(res$ances) %in% removed.strains
        if(any(temp)) {
            res$ances[temp] <<- NA
        }

        return(NULL)
    }



    ## MAIN SUB-FUNCTION: EXPANDING FROM ONE SEQUENCE - NON SPATIAL ##
    expand.one.strain <- function(seq, date, idx){
        toExpand[idx] <<- FALSE # this one is no longer to expand
        nbDes <- nb.desc()
        if(nbDes==0) return(NULL) # stop if no descendant
        newDates <- sapply(1:nbDes, function(i) date.dupli(date)) # find dates for descendants
        newDates <- newDates[newDates <= t.max] # don't store future sequences
        nbDes <- length(newDates)
        if(nbDes==0) return(NULL) # stop if no suitable date
        newSeq <- lapply(1:nbDes, function(i) seq.dupli(seq, newDates[i]-date)) # generate new sequences
        class(newSeq) <- "DNAbin" # lists of DNAbin vectors must also have class "DNAbin"
        newSeq <- as.matrix(newSeq) # list DNAbin -> matrix DNAbin with nbDes rows
        rownames(newSeq) <- seqname.gen(nbDes) # find new labels for these new sequences
        res$seq <<- rbind(res$seq, newSeq) # append to general output
        res$dates <<- c(res$dates, newDates) # append to general output
        res$ances <<- c(res$ances, rep(rownames(res$seq)[idx], nbDes)) # append to general output
        toExpand <<- c(toExpand, rep(TRUE, nbDes))
        return(NULL)
    }


    ## 2nd MAIN SUB-FUNCTION: EXPANDING FROM ONE SEQUENCE - SPATIAL ##
    expand.one.strain.xy <- function(seq, date, idx, cur.xy){
        toExpand[idx] <<- FALSE # this one is no longer to expand
        nbDes <- nb.desc()
        if(nbDes==0) return(NULL) # stop if no descendant
        newDates <- sapply(1:nbDes, function(i) date.dupli(date)) # find dates for descendants
        newDates <- newDates[newDates <= t.max] # don't store future sequences
        nbDes <- length(newDates)
        if(nbDes==0) return(NULL) # stop if no suitable date
        newSeq <- lapply(1:nbDes, function(i) seq.dupli(seq, newDates[i]-date)) # generate new sequences
        class(newSeq) <- "DNAbin" # lists of DNAbin vectors must also have class "DNAbin"
        newSeq <- as.matrix(newSeq) # list DNAbin -> matrix DNAbin with nbDes rows
        rownames(newSeq) <- seqname.gen(nbDes) # find new labels for these new sequences
        res$seq <<- rbind(res$seq, newSeq) # append to general output
        res$dates <<- c(res$dates, newDates) # append to general output
        res$ances <<- c(res$ances, rep(rownames(res$seq)[idx], nbDes)) # append to general output
        res$xy <<- rbind(res$xy, xy.dupli(cur.xy, nbDes))
        toExpand <<- c(toExpand, rep(TRUE, nbDes))
        return(NULL)
    }



    ## PERFORM SIMULATIONS - NON SPATIAL CASE ##
    if(!geo.sim){
        ## initialization
        res$seq <- matrix(rep(seq.gen(), ini.n), byrow=TRUE, nrow=ini.n)
        class(res$seq) <- "DNAbin"
        rownames(res$seq) <- 1:ini.n
        res$dates[1:ini.n] <- rep(0,ini.n)
        res$ances[1:ini.n] <- rep(NA,ini.n)
        toExpand <- rep(TRUE,ini.n)

        ## simulations: isn't simplicity beautiful?
        while(any(toExpand)){
            idx <- min(which(toExpand))
            expand.one.strain(res$seq[idx,], res$dates[idx], idx)
            resize.result()
        }


        ## SHAPE AND RETURN OUTPUT ##
        res$id <- as.character(1:length(res$ances))
        res$ances <- as.character(res$ances)
        names(res$dates) <- rownames(res$seq)
        res$call <- match.call()
        class(res) <- "haploGen"
        return(res)

    } # END NON-SPATIAL SIMULATIONS



    ## PERFORM SIMULATIONS - SPATIAL CASE ##
    if(geo.sim){
        ## some checks
        if(!is.null(mat.connect)) {
            if(nrow(mat.connect) != ncol(mat.connect)) stop("mat.connect is not a square matrix")
            if(nrow(mat.connect) != grid.size^2) stop("dimension of mat.connect does not match grid size")
        }

        ## initialization
        res$seq <- matrix(rep(seq.gen(), ini.n), byrow=TRUE, nrow=ini.n)
        class(res$seq) <- "DNAbin"
        rownames(res$seq) <- 1:ini.n
        res$dates[1:ini.n] <- rep(0,ini.n)
        res$ances[1:ini.n] <- rep(NA,ini.n)
        toExpand <- rep(TRUE,ini.n)

        if(is.null(ini.xy)){
            locStart <- xy.gen()
        } else{
            locStart <- as.vector(ini.xy)[1:2]
        }
        res$xy <- matrix(rep(locStart, ini.n), byrow=TRUE, nrow=ini.n)
        colnames(res$xy) <- c("x","y")

        ##cat("nb.strains","iteration.time",file="haploGenTime.out") # for debugging


        ## simulations: isn't simplicity beautiful?
        while(any(toExpand)){
            ##time.previous <- Sys.time() # FOR DEBUGGING
            idx <- min(which(toExpand))
            expand.one.strain.xy(res$seq[idx,], res$dates[idx], idx, res$xy[idx,])
            resize.result.xy()
            ## VERBOSE OUTPUT FOR DEBUGGING ##
            ## cat("\nNb strains:",length(res$ances),"/",max.nb.haplo)
            ##             cat("\nLatest date:", max(res$dates),"/",t.max)
            ##             cat("\nRemaining strains to duplicate", sum(toExpand))
            ##             cat("\n",append=TRUE,file="haploGenTime.out")
            ##             iter.time <- as.numeric(difftime(Sys.time(),time.previous,unit="sec"))
            ##             time.previous <- Sys.time()
            ##             cat(c(length(res$ances), iter.time),append=
            ##TRUE,file="haploGenTime.out")
        ## END DEBUGGING VERBOSE ##
        }

        ## VERBOSE OUTPUT FOR DEBUGGING ##
        ##cat("\nSimulation time stored in haploGenTime.out\n")

        ## SHAPE AND RETURN OUTPUT ##
        res$id <- as.character(1:length(res$ances))
        res$ances <- as.character(res$ances)
        names(res$dates) <- rownames(res$seq)

        class(res) <- "haploGen"
        res$call <- match.call()
        return(res)

    } # end SPATIAL SIMULATIONS


} # end haploGen








##################
## print.haploGen
##################
print.haploGen <- function(x, ...){

    cat("\t\n========================")
    cat("\t\n= simulated haplotypes =")
    cat("\t\n=  (haploGen object)   =")
    cat("\t\n========================\n")

    cat("\nSize :", length(x$ances),"haplotypes")
    cat("\nHaplotype length :", ncol(x$seq),"nucleotids")
    cat("\nProportion of NA ancestors :", signif(mean(is.na(x$ances)),5))
    cat("\nNumber of known ancestors :", sum(!is.na(x$ances)))
    nbAncInSamp <- sum(x$ances %in% labels(x))
    cat("\nNumber of ancestors within the sample :", nbAncInSamp)
    cat("\nDate range :", min(x$dates,na.rm=TRUE),"-",max(x$dates,na.rm=TRUE))
    ##nUniqSeq <- length(unique(apply(as.character(x$seq),1,paste,collapse="")))
    ##cat("\nNumber of unique haplotypes :", nUniqSeq)

    cat("\n\n= Content =")
    for(i in 1:length(x)){
        cat("\n")

        cat(paste("$", names(x)[i], sep=""),"\n")
        if(names(x)[i] %in% c("seq","call")) {
            print(x[[i]])
        } else if(names(x)[i]=="xy"){
            print(head(x[[i]]))
            if(nrow(x[[i]]>6)) cat("    ...\n")
        } else cat(head(x[[i]],6), ifelse(length(x[[i]])>6,"...",""),"\n")
    }


    return(NULL)
} # end print.haploGen






##############
## [.haploGen
##############
"[.haploGen" <- function(x,i,j,drop=FALSE){
    res <- x
    res$seq <- res$seq[i,,drop=FALSE]
    res$id <- res$id[i]
    res$ances <- res$ances[i]
    res$ances[!res$ances %in% res$id] <- NA
    res$dates <- res$dates[i]
    if(!is.null(res$xy)) res$xy <- res$xy[i,,drop=FALSE]

    return(res)
}






##################
## labels.haploGen
##################
labels.haploGen <- function(object, ...){
    return(object$id)
}





#######################
## as.POSIXct.haploGen
#######################
as.POSIXct.haploGen <- function(x, tz="", origin=as.POSIXct("2000/01/01"), ...){
    res <- as.POSIXct(x$dates*24*3600, origin=origin)
    return(res)
}






#####################
## seqTrack.haploGen
#####################
seqTrack.haploGen <- function(x, best=c("min","max"), prox.mat=NULL, ...){
    myX <- dist.dna(x$seq, model="raw")
    x.names <- labels(x)
    x.dates <- as.POSIXct(x)
    seq.length <- ncol(x$seq)
    myX <- myX * seq.length
    myX <- as.matrix(myX)
    prevCall <- as.list(x$call)
    if(is.null(prevCall$mu)){
        mu0 <- 0.0001
    } else {
        mu0 <- eval(prevCall$mu)
    }
    res <- seqTrack(myX, x.names=x.names, x.dates=x.dates, best=best, prox.mat=prox.mat,...)
    return(res)
}






########################
## as.seqTrack.haploGen
########################
as.seqTrack.haploGen <- function(x){
    ## x.ori <- x
    ## x <- na.omit(x)
    toSetToNA <- x$dates==min(x$dates)
    res <- list()
    res$id <- labels(x)
    res <- as.data.frame(res)
    res$ances <- x$ances
    res$ances[toSetToNA] <- NA
    res$weight <- 1 # ??? have to recompute that...
    res$weight[toSetToNA] <- NA
    res$date <- as.POSIXct(x)[labels(x)]
    res$ances.date <- as.POSIXct(x)[x$ances]

    ## set results as indices rather than labels
    res$ances <- match(res$ances, res$id)
    res$id <- 1:length(res$id)

    ## SET CLASS
    class(res) <- c("seqTrack", "data.frame")

    return(res)
}






################
## plotHaploGen
################
plotHaploGen <- function(x, annot=FALSE, date.range=NULL, col=NULL, bg="grey", add=FALSE, ...){

    ## SOME CHECKS ##
    if(class(x)!="haploGen") stop("x is not a haploGen object")
    if(is.null(x$xy)) stop("x does not contain xy coordinates - try converting to graphNEL for plotting")


    ## ## CONVERSION TO A SEQTRACK-LIKE OBJECT ##
    xy <- x$xy
    res <- as.seqTrack.haploGen(x)

    ##     res <- list()
    ##     res$id <- labels(x)
    ##     res <- as.data.frame(res)
    ##     res$ances <- x$ances
    ##     res$ances[toSetToNA] <- NA
    ##     res$weight <- 1 # ??? have to recompute that...
    ##     res$weight[toSetToNA] <- NA
    ##     res$date <- as.POSIXct(x.ori)[labels(x)]
    ##     res$ances.date <- as.POSIXct(x.ori)[x$ances]
    ##     ## set results as indices rather than labels
    ##     res$ances <- match(res$ances, res$id)
    ##     res$id <- 1:length(res$id)


    ## CALL TO PLOTSEQTRACK ##
    plotSeqTrack(x=res, xy=xy, annot=annot, date.range=date.range,
                        col=col, bg=bg, add=add, ...)

    return(invisible(res))

} # end plotHaploGen






###################
## sample.haploGen
###################
sample.haploGen <- function(x, n){
##sample.haploGen <- function(x, n, rDate=.rTimeSeq, arg.rDate=NULL){
    ## EXTRACT THE SAMPLE ##
    res <- x[sample(1:nrow(x$seq), n, replace=FALSE)]


    ## RETRIEVE SOME PARAMETERS FROM HAPLOSIM CALL
    prevCall <- as.list(x$call)
    if(!is.null(prevCall$mu)){
        mu0 <- eval(prevCall$mu)
    } else {
        mu0 <- 1e-04
    }

    if(!is.null(prevCall$seq.length)){
        L <- eval(prevCall$seq.length)
    } else {
        L <- 1000
    }

    ## truedates <- res$dates
    ## daterange <- diff(range(res$dates,na.rm=TRUE))

    ## if(identical(rDate,.rTimeSeq)){
    ##     sampdates <- .rTimeSeq(n=length(truedates), mu=mu0, L=L, maxNbDays=daterange/2)
    ## } else{
    ##     arg.rDate$n <- n
    ##     sampdates <- do.call(rDate, arg.rDate)
    ## }
    ## sampdates <- truedates + abs(sampdates)

    return(res)
} # end sample.haploGen






######################
## as.igraph.haploGen
######################
as.igraph.haploGen <- function(x, col.pal=redpal, ...){
    ## if(!require(igraph)) stop("package igraph is required for this operation")
    ## if(!require(ape)) stop("package ape is required for this operation")

    ## GET DAG ##
    from <- x$ances
    to <- x$id
    isNotNA <- !is.na(from) & !is.na(to)
    dat <- data.frame(from,to,stringsAsFactors=FALSE)[isNotNA,,drop=FALSE]
    vnames <- as.character(unique(unlist(dat)))
    out <- graph.data.frame(dat, directed=TRUE, vertices=data.frame(names=vnames, dates=x$dates[vnames]))

    ## SET WEIGHTS ##
    D <- as.matrix(dist.dna(x$seq,model="raw")*ncol(x$seq))
    temp <- mapply(function(i,j) return(D[i,j]), as.integer(from), as.integer(to))
    E(out)$weight <- temp[isNotNA]


    ## DATES FOR VERTICES
    V(out)$dates <- x$date

    ## SET EDGE LABELS ##
    E(out)$label <- E(out)$weight

    ## SET EDGE COLORS
    E(out)$color <- num2col(E(out)$weight, col.pal=col.pal, reverse=TRUE)

    ## SET LAYOUT ##
    ypos <- V(out)$dates
    ypos <- abs(ypos-max(ypos))
    attr(out, "layout") <- layout.fruchterman.reingold(out, params=list(miny=ypos, maxy=ypos))

    return(out)
} # end as.igraph.haploGen






#################
## plot.haploGen
#################
plot.haploGen <- function(x, y=NULL, col.pal=redpal, ...){
    ## if(!require(igraph)) stop("igraph is required")

    ## get graph ##
    g <- as.igraph(x, col.pal=col.pal)

    ## make plot ##
    plot(g, layout=attr(g,"layout"), ...)

    ## return graph invisibly ##
    return(invisible(g))

} # end plot.haploGen
















##########################
## as("haploGen", "graphNEL")
##########################
## if(require(graph)){
##     setOldClass("haploGen")

##     setAs("haploGen", "graphNEL", def=function(from){
##         if(!require(ape)) stop("package ape is required")
##         if(!require(graph)) stop("package graph is required")

##         N <- length(from$ances)
##         areNA <- is.na(from$ances)

##         ## EXTRACT WEIGHTS (nb of mutations)
##         M <- as.matrix(dist.dna(from$seq, model="raw")*ncol(from$seq))
##         rownames(M) <- colnames(M) <- from$id
##         w <- mapply(function(i,j) {M[i, j]}, i=from$ances[!areNA], j=from$id[!areNA])


##         ## CONVERT TO GRAPH
##         res <- ftM2graphNEL(ft=cbind(from$ances[!areNA], from$id[!areNA]), W=w, edgemode = "directed", V=from$id)
##         return(res)
##     })

## }






















## #####################
## ## seqTrackG.haploGen
## #####################
## seqTrackG.haploGen <- function(x, optim=c("min","max"), ...){
##     myX <- dist.dna(x$seq, model="raw")
##     x.names <- labels(x)
##     x.dates <- as.POSIXct(x)
##     seq.length <- ncol(x$seq)
##     myX <- myX * seq.length
##     prevCall <- as.list(x$call)
##     if(is.null(prevCall$mu)){
##         mu0 <- 0.0001
##     } else {
##         mu0 <- eval(prevCall$mu)
##     }
##     res <- seqTrackG(myX, x.names=x.names, x.dates=x.dates, best=optim,...)
##     return(res)
## }






##############################
## optimize.seqTrack.haploGen
##############################
## optimize.seqTrack.haploGen <- function(x, thres=0.2, optim=c("min","max"),
##                                        typed.chr=NULL, mu0=NULL, chr.length=NULL,
##                                        prox.mat=NULL, nstep=10, step.size=1e3,
##                                        rDate=.rTimeSeq, arg.rDate=NULL, rMissDate=.rUnifTimeSeq, ...){

##     x.names <- labels(x)
##     x.dates <- as.POSIXct(x)
##     seq.length <- ncol(x$seq)
##     myX <- dist.dna(x$seq, model="raw") * seq.length
##     prevCall <- as.list(x$call)
##     if(is.null(prevCall$mu)){
##         mu0 <- 0.0001
##     } else {
##         mu0 <- eval(prevCall$mu)
##     }

##     res <- optimize.seqTrack.default(x=myX, x.names=x.names, x.dates=x.dates,
##                                      typed.chr=typed.chr, mu0=mu0, chr.length=chr.length,
##                                      thres=thres, optim=optim, prox.mat=prox.mat,
##                                      nstep=nstep, step.size=step.size,
##                                      rDate=rDate, arg.rDate=arg.rDate, rMissDate=rMissDate, ...)
## } # end optimize.seqTrack.haploGen
