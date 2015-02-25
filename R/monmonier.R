# Algorithm to detect boundaries, based on Monmonier's algorithm
# Extended to any connection network
# Thibaut Jombart 2006-2008 (t.jombart@imperial.ac.uk)



#####################
# function monmonier
#####################
monmonier <- function(xy, dist, cn, threshold=NULL, bd.length=NULL, nrun=1,
                      skip.local.diff=rep(0,nrun), scanthres=is.null(threshold), allowLoop=TRUE){
if(!require(spdep)) stop("The package spdep is required but not installed")
if(!inherits(cn,"nb")) stop('cn is not a nb object')
if(is.data.frame(xy)) xy <- as.matrix(xy)
if(!is.matrix(xy)) stop('xy must be a matrix')
if(!inherits(dist,"dist")) stop('Argument \'dist\' must be a distance matrix of class dist')
if(nrow(xy) != nrow(as.matrix(dist))) stop('Number of sites and number of observations differ')

## set to TRUE to debug
## DEBUG <- FALSE

## if(DEBUG) {
##     plot.nb(cn, xy, col="grey", points=FALSE)
##     if(exists("x1")) {
##         x1 <<- x1
##         s.value(xy,x1,add.p=TRUE)
##     }
##     text(xy,lab=1:nrow(xy),font=2,col="darkgreen")
## }

## PRECISION of the xy coordinates (in digits)
## used when coordinates are inputed in C code.
PRECISION=12

## conversion of the connection network
cn.nb <- cn
cn <- nb2neig(cn)
## binary matrix of neighbourhood
M <- neig2mat(cn)
## distance matrix
D <- as.matrix(dist)
## matrix of distances among neighbours
D <- M*D

## set/check the threshold value of the distances among neighbours
## default: the median of all distances among neighbours
if(is.null(threshold) || threshold<0) {Dlim <- summary(unique(D[D>0]))[5]} else {Dlim <- threshold}

## handle 'bd.length' / must prevail over threshold.
if(!is.null(bd.length)) {
    if(bd.length < 1) bd.length <- 1
    Dlim <- 0
    scanthres <- FALSE
} else {
    bd.length <- nrow(xy)^2
}


if(scanthres){
  plot(sort(unique(D)[unique(D) > 0],decreasing=TRUE),main="Local distances plot",type="l",xlab="rank",ylab="Sorted local distances")
  abline(h=Dlim,lty=2)
  mtext("Dashed line indicates present threshold")
  cat("Indicate the threshold (\'d\' for default): ")
  temp <- as.character(readLines(n = 1))
  if(toupper(temp)!="D") { Dlim <- as.numeric(temp) }
}

#### OLD VERSION ####
## ## build a data.frame of coords of neighbours, with columns x1 y1 x2 y2
## listCpl <- neig.util.GtoL(neig2mat(cn))
## allSeg <- cbind(xy[,1][listCpl[,1]] , xy[,2][listCpl[,1]] , xy[,1][listCpl[,2]] , xy[,2][listCpl[,2]])
## colnames(allSeg) <- c('xP','yP','xQ','yQ')
#### END OLD VERSION ####

## Compute once and for all a matrix matSegVal with columns A, B and val
## A: index of one vertex
## B: index if another vertex, connected (beighbour) to B
## val: distance between A and B
## This table is such that "val" is sorted in decreasing order.
matSegVal <- which(D>0,arr.ind=TRUE)
matSegVal.names <- apply(matSegVal,1,function(vec) paste(sort(vec),collapse="-"))
rownames(matSegVal) <- matSegVal.names # used to get unique segments (cause D is symmetric)
matSegVal <- matSegVal[unique(rownames(matSegVal)),]
temp <- apply(matSegVal,1,function(vec) D[vec[1],vec[2]])
matSegVal <- cbind(matSegVal,temp) # matSegVal has its 3 colums
idx <- order(matSegVal[,3], decreasing=TRUE)
matSegVal <- matSegVal[idx,]

## allSeg is a matrix giving the coordinates of all edges (in row), x1 y1 x2 y2.
allSeg <- cbind(xy[matSegVal[,1],] , xy[matSegVal[,2],])
rownames(allSeg) <- rownames(matSegVal)
colnames(allSeg) <- c("x1","y1","x2","y2")

## auxiliary function used to remove a segment from matSegVal (A,B: indices of vertices of the removed edge)
rmFromMatSegVal <- function(A,B){
    AB.name <- paste(sort(c(A,B)), collapse="-")
    toRemove <- match(AB.name, rownames(matSegVal))
    matSegVal <<- matSegVal[-toRemove,]
}

## Given a rank 'rang', getNext retrieves the two corresponding neighbours (A and B),
## their middle (M), and the corresponding value.
## A is the index of the vertex, xA and yA are its coordinates
getNext <- function(rang){
    ## ## ## OLD VERSION ## ## ##
    ## val <- unique(sort(D,decreasing=TRUE))[rang]
    ##  A <- which(round(D,10)==round(val,10),TRUE)[1,1] # beware: == used on doubles
    ##     B <- which(round(D,10)==round(val,10),TRUE)[1,2] # beware: == used on doubles
    ##     xA <- xy[A,1]
    ##     yA <- xy[A,2]
    ##     xB <- xy[B,1]
    ##     yB <- xy[B,2]
    ## 	xM <- (xA + xB)/2
    ##     yM <- (yA + yB)/2
    ## ## ## END OLD VERSION ## ## ##

    #### Must take into account that the some local distances occur for different segments
    #### Create a table with columns A, B and val

    temp <- matSegVal[rang,]
    A <- temp[1]
    B <- temp[2]
    val <- temp[3]
    xA <- xy[A,1]
    yA <- xy[A,2]
    xB <- xy[B,1]
    yB <- xy[B,2]
    xM <- mean(c(xA,xB))
    yM <- mean(c(yA,yB))

    return( list(A=c(A,xA,yA), B=c(B,xB,yB), M=c(xM,yM), val=val) )
}

## checkNext returns TRUE if it is ok to draw the segment MN, i.e. if MN does not cross another edge,
## and FALSE otherwise
## - M and N define the segment of interest
## M = c(xM,yM) ; N=c(xN,yN)
## -segMat is the matrix of all known edges
## A,B,C,D are not coordinates, but indices of vertices
## -AB is the edge whose middle is M ; given to avoid checking MN vs AB
## -CD is the edge whose middle is N ; given to avoid checking MN vs CD
## - toRemove: a vector of characters naming other edges in segMat to be removed before checking
checkNext <- function(M,N,A,B,C,D,segMat=allSeg,toRemove=NULL,curDir=NULL){
    ## orientation of the segment MN
    xmin <- min(M[1],N[1])
    xmax <- max(M[1],N[1])
    ymin <- min(M[2],N[2])
    ymax <- max(M[2],N[2])

### From here segment that would be difficult to check for crossing are removed
### If MN is the segment of interest, with M middle of AB and N middle of CD,

    subsetSeg <- segMat
    ## The segment that was just drawn before is removed from the checks for crossing (messy code 2)
    if(!is.null(curDir)){
        prevSeg <- grep(paste("dir",curDir,sep=""),rownames(subsetSeg))
        if(length(prevSeg)>0){ # this is not the 1st seg. in this dir
            prevSeg <- prevSeg[length(prevSeg)]
            subsetSeg <- segMat[-prevSeg,,drop=FALSE]
        } else { # this is the 1st seg. in this direction -> rm 1st seg of other dir
            otherDir <- ifelse(curDir==1L, 2, 1)
            prevSeg <- grep(paste("dir",otherDir,"-1-2",sep=""),rownames(subsetSeg))
            if(length(prevSeg)>0) { subsetSeg <- segMat[-prevSeg,,drop=FALSE] }
        }
    }

    ## AB and CD are now taken out (doubtful code 2 returned by the C procedure CheckAllSeg)
    AB.name <- paste(sort(c(A,B)),collapse="-")
    AB.idx <- match(AB.name, rownames(subsetSeg))
    CD.name <- paste(sort(c(C,D)),collapse="-")
    CD.idx <- match(CD.name, rownames(subsetSeg))
    subsetSeg <- subsetSeg[-c(AB.idx,CD.idx),,drop=FALSE]

    ## subsetSeg is a matrix, each of its rows being a segment: xP,yP,xQ,yQ
    ## edges out of the square of diagonal MN cannot be crossed, so they are removed from
    ## the checks for crossing
    subsetSeg <- subsetSeg[!(subsetSeg[,1] < xmin & subsetSeg[,3] < xmin),,drop=FALSE]
    subsetSeg <- subsetSeg[!(subsetSeg[,1] > xmax & subsetSeg[,3] > xmax),,drop=FALSE]
    subsetSeg <- subsetSeg[!(subsetSeg[,2] < ymin & subsetSeg[,4] < ymin),,drop=FALSE]
    subsetSeg <- subsetSeg[!(subsetSeg[,2] > ymax & subsetSeg[,4] > ymax),,drop=FALSE]

    ## handle toRemove here
    if(!is.null(toRemove)){
        idx <- match(toRemove, rownames(subsetSeg))
        idx <- idx[!is.na(idx)]
        if(length(idx)>0) { subsetSeg <- subsetSeg[-idx,,drop=FALSE] }
    }

    ## temp is used for the output of CheckAllSeg
    ## initialized at 10, which is never returned by CheckAllSeg
    temp <- as.integer(10)

    ## call to CheckAllSeg
    ## =======================
    ## output code for CheckAllSeg:
    ## - 0: no intersection
    ## - 1: all kind of intersection, including the codes from SegSeg function
    ## (a C function inside monmonier-utils.C):
    ## - 3 : The segments collinearly overlap, sharing at least a point.
    ## - 2 : An endpoint (vertex) of one segment is on the other segment,
    ## but segments aren't collinear.
    ## - 1 : The segments intersect properly (i.e. not case 2 or 3)
    ## =======================
    ##
    ## restore the matrix type if there is only one segment to check for crossing in subsetSeg

    ## round down coordinates in subsetSeg
    subsetSeg <- round(subsetSeg, digits=PRECISION)

    if(nrow(subsetSeg)>0) {
        temp <- .C("CheckAllSeg",as.integer(nrow(subsetSeg)),as.integer(ncol(subsetSeg)),
                   as.double(as.matrix(subsetSeg)), as.double(M), as.double(N), temp,PACKAGE="adegenet")[[6]]
    } else {temp <- 0}

    ## for debugging
    ##  if(DEBUG) {
    ##                 if(temp==1)  cat("\n can't go there (code",temp,")") else cat("\n new segment ok (code",temp,")")
    ##             }

    ## if a code 1 or 3 is returned, CheckAllSeg returns FALSE
    ## else it returns TRUE
    ## additional control used (code 10, CheckAllSeg failed)
    if(temp==10) stop("CheckAllSeg failure (returned value=10, i.e. unchanged, not computed). Please report the bug.")
    if(temp==1 | temp==2 | temp==3) return(FALSE) else return(TRUE)
} # end of checkNext

## result object is created
## this is a list with one component per run (a run = a boundary)
result <-list()
for(run in 1:nrun){
    result[[run]] <- list(dir1=list(),dir2=list())
}


#### MAIN FUNCTION HERE
#### a for loop is used to handle several runs
## each boundary seeks one starting point, and then expands both sides (dir1 / dir2)
for(run in 1:nrun){

    ## handle skip.local.diff here
    ## the corresponding values of distance are set definitely to -1
    ##   if(skip.local.diff[run] >0)
    ##         for(i in 1:skip.local.diff[run]){
    ##             temp <- getNext(1)
    ##             ## D[temp$A[1],temp$B[1] ] <- -1
    ##             ## D[temp$B[1],temp$A[1] ] <- -1
    ##             rmFromMatSegVal(temp$A[1],temp$B[1])
    ##         }

    ## starting point: find the highest distance among neighbours
    ## then expand
    currentDir1 <- getNext(1 + skip.local.diff[run])
    currentDir2 <- currentDir1
    current.bd.length <- 1
    if(currentDir1$val<=Dlim) stop(paste('Algorithm reached the threshold value at the first step of run',run))
    result[[run]]$dir1[[1]] <- currentDir1 # first point dir1
    result[[run]]$dir2[[1]] <- currentDir2 # first point dir2 (same as dir1)
    ## D[result[[run]]$dir1[[1]]$A[1],result[[run]]$dir1[[1]]$B[1]] <- -1 # update D matrix
    ## D[result[[run]]$dir1[[1]]$B[1],result[[run]]$dir1[[1]]$A[1]] <- -1 # update D matrix
    rmFromMatSegVal(result[[run]]$dir1[[1]]$A[1], result[[run]]$dir1[[1]]$B[1]) # update matSegVal

    ## dir1 => i1: result index; s1: index for the rank of the distance between neighbours (decreasing order)
    ## dir2 => i2: result index; s2: index for the rank of the distance between neighbours (decreasing order)
    i1 <- 1
    s1 <- 1
    i2 <- 1
    s2 <- 2
    n <- nrow(D)

    ## logical handling the expansion of a boundary
    keepExpanding <- ((current.bd.length < bd.length)
                      && ((currentDir1$val>Dlim)|(currentDir2$val>Dlim))
                      && (s1 < nrow(matSegVal) | s2 < nrow(matSegVal)) )

    #### This while loop has the following behavior:
    ## as long as the keepExpanding is true, we try to expand the boundary by
    ## 1) finding the highest next distance among neighbours
    ## 2) test if we can draw the corresponding segment
    ## 3) - if we can, store the result, increment i1 or i2, reset s1 and s2, erase the edge
    ## 3) - if we can't, take the next distance among neighbours (incrementing s1 and s2)
    ## 4) get back to 1)
    while(keepExpanding){
        hasExpanded <- FALSE # used to test if it is relevant to check for looping
        ##  if(DEBUG){
        ##                             points(currentDir1$M[1],currentDir1$M[2],col="white",pch=20)
        ##                             points(currentDir2$M[1],currentDir2$M[2],col="white",pch=20)
        ##                         }
        if(s1 <= nrow(matSegVal)) { currentDir1 <- getNext(s1) }
        if(s2 <= nrow(matSegVal)) { currentDir2 <- getNext(s2) }
        ##  if(DEBUG){
        ##                             cat("\n\n ## dir1: trying edge",currentDir1$A[1],"-",currentDir1$B[1])
        ##                             points(currentDir1$M[1],currentDir1$M[2],col="blue",pch=20)
        ##                             readline("\npress enter")
        ##                         }

        ## first direction (dir1)
        if( currentDir1$val > Dlim && s1 <= nrow(matSegVal)) {
            if(checkNext(result[[run]]$dir1[[length(result[[run]]$dir1)]]$M,
                         currentDir1$M,
                         result[[run]]$dir1[[length(result[[run]]$dir1)]]$A[1],
                         result[[run]]$dir1[[length(result[[run]]$dir1)]]$B[1],
                         currentDir1$A[1],
                         currentDir1$B[1], curDir=as.integer(1))) {
                i1 <- i1+1
                result[[run]]$dir1[[i1]] <- currentDir1
                ## update the matrix of diff. between neighbours
                ## D[result[[run]]$dir1[[i1]]$A[1],result[[run]]$dir1[[i1]]$B[1]] <- -1
                ## D[result[[run]]$dir1[[i1]]$B[1],result[[run]]$dir1[[i1]]$A[1]] <- -1
                rmFromMatSegVal(result[[run]]$dir1[[i1]]$A[1],result[[run]]$dir1[[i1]]$B[1] ) # update matSegVal

                s1 <- 1
                ## update existing segments
                allSeg <- rbind(allSeg,c(result[[run]]$dir1[[i1-1]]$M,result[[run]]$dir1[[i1]]$M) )
                rownames(allSeg) <- c(rownames(allSeg)[-nrow(allSeg)] , paste("dir1",i1-1,i1,sep="-"))
                ## add 1 to the boundary length
                current.bd.length <- current.bd.length + 1
                hasExpanded <- TRUE
               ##  if(DEBUG) {
                ##                                                     arrows(result[[run]]$dir1[[i1-1]]$M[1], result[[run]]$dir1[[i1-1]]$M[2],
                ##                                                            result[[run]]$dir1[[i1]]$M[1], result[[run]]$dir1[[i1]]$M[2],col="blue")
                ##                                                 }

            } else{ s1 <- s1+1 }
        } # end "if( currentDir1$val>Dlim)"

          ## if(DEBUG){
        ##                             cat("\n\n ## dir2: trying edge",currentDir2$A[1],"-",currentDir2$B[1])
        ##                             points(currentDir2$M[1],currentDir2$M[2],col="red",pch=20)
        ##                             readline("\npress enter")
        ##                         }

        ## second direction (dir2)
        if( currentDir2$val > Dlim  && s2 <= nrow(matSegVal)) {
            if(checkNext(result[[run]]$dir2[[length(result[[run]]$dir2)]]$M,
                         currentDir2$M,
                         result[[run]]$dir2[[length(result[[run]]$dir2)]]$A[1],
                         result[[run]]$dir2[[length(result[[run]]$dir2)]]$B[1],
                         currentDir2$A[1],
                         currentDir2$B[1], curDir=as.integer(2))) {
                i2 <- i2+1
                result[[run]]$dir2[[i2]] <- currentDir2
                ## update the matrix of diff. between neighbours
                ## D[result[[run]]$dir2[[i2]]$A[1],result[[run]]$dir2[[i2]]$B[1]] <- -1
                ## D[result[[run]]$dir2[[i2]]$B[1],result[[run]]$dir2[[i2]]$A[1]] <- -1
                rmFromMatSegVal(result[[run]]$dir2[[i2]]$A[1],result[[run]]$dir2[[i2]]$B[1])
                s2 <- 1
                ## update existing segments
                allSeg <- rbind(allSeg,c(result[[run]]$dir2[[i2-1]]$M,result[[run]]$dir2[[i2]]$M) )
                rownames(allSeg) <- c(rownames(allSeg)[-nrow(allSeg)] , paste("dir2",i2-1,i2,sep="-"))
                ## add 1 to the boundary length
                current.bd.length <- current.bd.length + 1
                hasExpanded <- TRUE
                 ## if(DEBUG){
                ##                                                     arrows(result[[run]]$dir2[[i2-1]]$M[1], result[[run]]$dir2[[i2-1]]$M[2],
                ##                                                            result[[run]]$dir2[[i2]]$M[1], result[[run]]$dir2[[i2]]$M[2],col="red",cex=2)
                ##                                                 }

            } else{ s2 <- s2+1 }
        } # end "if( currentDir2$val>Dlim)"

        ## update the logical for the while loop
        keepExpanding <- ((current.bd.length < bd.length)
                          && ((currentDir1$val>Dlim)|(currentDir2$val>Dlim))
                          && (s1 <= nrow(matSegVal) | s2 <= nrow(matSegVal)) )

        ## handle the looping of a boundary
        if(hasExpanded && (current.bd.length>3) && allowLoop){
            ## check if the two ends of the boundary can be joined
            ## remove segments ending each direction (to avoid messy code 2 in checkNext)
            terminalEdges <- c(paste("dir1",i1-1,i1,sep="-"),
                               paste("dir2",i2-1,i2,sep="-"))

            canLoop <- checkNext(result[[run]]$dir1[[length(result[[run]]$dir1)]]$M,
                                 result[[run]]$dir2[[length(result[[run]]$dir2)]]$M,
                                 result[[run]]$dir1[[length(result[[run]]$dir1)]]$A[1],
                                 result[[run]]$dir1[[length(result[[run]]$dir1)]]$B[1],
                                 result[[run]]$dir2[[length(result[[run]]$dir2)]]$A[1],
                                 result[[run]]$dir2[[length(result[[run]]$dir2)]]$B[1],
                                 toRemove=terminalEdges
                                 )

            if(canLoop) {
                ## add the last, closing segment
                result[[run]]$dir1[[length(result[[run]]$dir1)+1]] <- result[[run]]$dir2[[length(result[[run]]$dir2)]]
                ## stop expanding
                keepExpanding <- FALSE
                ## update existing segments
                allSeg <- rbind(allSeg,c(result[[run]]$dir1[[length(result[[run]]$dir1)]]$M,
                                         result[[run]]$dir2[[length(result[[run]]$dir2)]]$M) )
            } # end looping of the boundary

        }

        ## output for debugging
        ##  if(DEBUG) {
        ##                     cat("\n","s1:",s1,"s2:",s2,"i1:",i1,"i2:",i2,"D1:",
        ##                         currentDir1$val,"D2:",currentDir2$val,"Dlim:",Dlim,
        ##                         "nrow(matSegVal)",nrow(matSegVal),"\n")
        ##                     cat("\n","D1:",currentDir1$val,"D2:",currentDir2$val,"Dlim:",Dlim,
        ##                         "cur.bd.le:", current.bd.length,"max length:", bd.length,"\n",
        ##                         "s1:",s1,"s2:",s2,"maxS:",nrow(matSegVal))
        ##                 }

    } # end of one given run

} # end for all run

# build the final output
# this is a list of class monmonier
# each element correspond to a run, i.e. to a potential boundary
# the output also contains the number of runs ($nrun) and the matched call ($call).
output=list()
for(run in 1:nrun){
	runname <- paste('run',run,sep='')
	output[[runname]] <- list(dir1=list(),dir2=list())
	# dir 1 #
	output[[runname]]$dir1$path <- matrix(-1, ncol=2,nrow=length(result[[run]]$dir1))
	colnames(output[[runname]]$dir1$path) <- c('x','y')
	rownames(output[[runname]]$dir1$path) <- paste('Point',1:nrow(output[[run]]$dir1$path),sep='_')

	for(i in 1:length(result[[run]]$dir1)) {
		output[[runname]]$dir1$path[i,] <- result[[run]]$dir1[[i]]$M
		output[[runname]]$dir1$values[i] <- result[[run]]$dir1[[i]]$val
	}

        # dir 2 #
	output[[runname]]$dir2$path <- matrix(-1, ncol=2,nrow=length(result[[run]]$dir2))
	colnames(output[[runname]]$dir2$path) <- c('x','y')
	rownames(output[[runname]]$dir2$path) <- paste('Point',1:nrow(output[[run]]$dir2$path),sep='_')
        for(i in 1:length(result[[run]]$dir2)) {
		output[[runname]]$dir2$path[i,] <- result[[run]]$dir2[[i]]$M
		output[[runname]]$dir2$values[i] <- result[[run]]$dir2[[i]]$val
	}

}

output$nrun <- nrun
output$threshold <- Dlim
output$xy <- xy
output$cn <- cn.nb
output$call <- match.call()
class(output) <- 'monmonier'
return(output)
}




##########################
# function plot.monmonier
##########################
plot.monmonier <- function(x, variable=NULL,displayed.runs=1:x$nrun,
                           add.arrows=TRUE, col='blue',lty=1,bwd=4, clegend=1,csize=0.7,
                           method = c('squaresize','greylevel'),sub='',csub=1,possub='topleft',
                           cneig=1,pixmap=NULL,contour=NULL,area=NULL,add.plot=FALSE,...){

    if (!inherits(x, "monmonier")) stop("Use only with 'monmonier' objects")
    if(!is.null(variable) & !is.numeric(variable)) stop('If provided, variable must be numeric.\n')


    xy <- x$xy
    cpoint <- 0

    if(cneig>0) {neig <- nb2neig(x$cn)} else {neig <- NULL}
    if(is.null(variable)){
	variable <- rep(1,nrow(xy))
	csize <- 0
	cpoint <- 1
	clegend <- 0
    }
    s.value(xy,variable,grid=FALSE,include.origin=FALSE,addaxes=FALSE,neig=neig,
            cneig=cneig,clegend=clegend,csize=csize,cpoint=cpoint,pch=20,pixmap=pixmap,
            method=method,sub=sub,csub=csub,possub=possub,add.plot=add.plot)
    opar <- par(no.readonly=TRUE)
    on.exit(par(mar=opar$mar))
    par(mar=c(0,0,0,0))

    for(run in displayed.runs){
        obj <- x[[run]]
        if(length(col)!=x$nrun) col <- rep(col,x$nrun)
        if(length(lty)!=x$nrun) lty <- rep(lty,x$nrun)
        if(length(obj$dir1$values) == 0) stop(paste('Monmonier object of run', run, 'is empty (no point in the path)\n'))

        if(length(obj$dir1$values) == 1 && length(obj$dir2$values) == 1) {
            points(obj$dir1$path[1],obj$dir1$path[2],pch=20,col=col[run],...)
        } else{
            ## handle boundary width
            ## the largest part corresponds to the highest distance among neighbours
            val.1 <- obj$dir1$values
            val.2 <- obj$dir2$values
            n1 <- length(val.1)
            n2 <- length(val.2)
            cex.bwd.1 <- ( val.1[1:(n1-1)] + val.1[2:n1] )/2
            cex.bwd.2 <- ( val.2[1:(n2-1)] + val.2[2:n2] )/2
            cex.bwd.max <- max(c(cex.bwd.1,cex.bwd.2), na.rm=TRUE)
            cex.bwd.1 <- cex.bwd.1/max(cex.bwd.max)
            cex.bwd.2 <- cex.bwd.2/max(cex.bwd.max)
            ## amplify the differences
            ## cex.bwd.1 <- cex.bwd.1^1.5
            ## cex.bwd.2 <- cex.bwd.2^1.5


            if(add.arrows) {
                if(n1>1) arrows(obj$dir1$path[1:(nrow(obj$dir1$path)-1),1],
                                obj$dir1$path[1:(nrow(obj$dir1$path)-1),2],
                                obj$dir1$path[2:nrow(obj$dir1$path),1],
                                obj$dir1$path[2:nrow(obj$dir1$path),2],
                                lwd=bwd*cex.bwd.1,angle=20,length=0.2,col=col[run],lty=lty[run],...)

                if(n2>1) arrows(obj$dir2$path[1:(nrow(obj$dir2$path)-1),1],
                                obj$dir2$path[1:(nrow(obj$dir2$path)-1),2],
                                obj$dir2$path[2:nrow(obj$dir2$path),1],
                                obj$dir2$path[2:nrow(obj$dir2$path),2],
                                lwd=bwd*cex.bwd.2,angle=20,length=0.2,col=col[run],lty=lty[run],...)
            } else {
                if(n1>1) segments(obj$dir1$path[1:(nrow(obj$dir1$path)-1),1],
                                  obj$dir1$path[1:(nrow(obj$dir1$path)-1),2],
                                  obj$dir1$path[2:nrow(obj$dir1$path),1],
                                  obj$dir1$path[2:nrow(obj$dir1$path),2],
                                  lwd=bwd*cex.bwd.1,col=col[run],lty=lty[run],...)

                if(n2>1)segments(obj$dir2$path[1:(nrow(obj$dir2$path)-1),1],
                                 obj$dir2$path[1:(nrow(obj$dir2$path)-1),2],
                                 obj$dir2$path[2:nrow(obj$dir2$path),1],
                                 obj$dir2$path[2:nrow(obj$dir2$path),2],
                                 lwd=bwd*cex.bwd.2,col=col[run],lty=lty[run],...)
            }
        } # end else
    } # end for

} # end function



#####################
# print function
#####################
print.monmonier <- function(x, ...){
cat("\t\n###########################################################")
cat("\t\n# List of paths of maximum differences between neighbours #")
cat("\t\n#           Using a Monmonier based algorithm             #")
cat("\t\n###########################################################\n")
cat('\n$call:')
print(x$call)

cat('\n      # Object content #')
cat("\nClass: ", class(x))
cat('\n$nrun (number of successive runs): ', x$nrun)
if(x$nrun==1)
cat('\n$run1: run of the algorithm')
else if(x$nrun==2)
cat('\n$run1, $run2: runs of the algorithm')
else cat('\n$run1 ... $run',x$nrun, ': runs of the algorithm',sep='')
cat('\n$threshold (minimum difference between neighbours): ', x$threshold)
cat("\n$xy: spatial coordinates")
cat("\n$cn: connection network")

cat('\n\n      # Runs content #')
for(i in 1:x$nrun){
	cat('\n# Run',i)
        # dir 1 #
        cat('\n# First direction')
        cat('\nClass: ', class(x$run1$dir1))
	cat('\n$path:\n')
	print(head(x[[i]]$dir1$path,n=3))
	if(nrow(x[[i]]$dir1$path) >3) cat('...\n')
	cat('\n$values:\n',head(x[[i]]$dir1$values,n=3))
	if(length(x[[i]]$dir1$values)>3) cat(' ...')
        # dir 2 #
        cat('\n# Second direction')
        cat('\nClass: ', class(x$run1$dir2))
	cat('\n$path:\n')
	print(head(x[[i]]$dir2$path,n=3))
	if(nrow(x[[i]]$dir2$path) >3) cat('...\n')
	cat('\n$values:\n',head(x[[i]]$dir2$values,n=3))
	if(length(x[[i]]$dir2$values)>3) cat(' ...')

	cat('\n')
	lenTheo <- x$nrun + 5
	if(length(names(x))> lenTheo) {
		cat('Other elements: \n')
		cat(names(x)[(lenTheo+1) : length(x)])
		}
	cat('\n')
	}
}



##############################
# function optimize.monmonier
##############################
optimize.monmonier <- function(xy,dist,cn,ntry=10, bd.length=NULL, return.best=TRUE,
                               display.graph=TRUE,threshold=NULL,scanthres=is.null(threshold),allowLoop=TRUE){

## move X to the arguments if we want to optimize an already created object
X <- NULL
#if( any(is.null(xy), is.null(dist),  is.null(cn)) & is.null(X) ) stop("Please provide either xy, dist and cn or a monmonier object (X)")

## if X is a monmonier object...
if(inherits(X,what="monmonier")){
  obj <- as.list(X$call)
  xy <- obj$xy
  dist <- obj$dist
  cn <- obj$cn
}

cn.nb <- cn
cn <- nb2neig(cn)
M <- neig2mat(cn)
D <- as.matrix(dist)
D <- M*D
if(is.null(threshold) || (threshold<0)) {Dlim <- summary(unique(D[D>0]))[5]} else {Dlim <- threshold}

if(scanthres){
  plot(sort(unique(D)[unique(D) > 0],decreasing=TRUE),main="Local distances plot", type="l", xlab="rank",ylab="Sorted local distances")
  abline(h=Dlim,lty=2)
  mtext("Dashed line indicates present threshold")
  cat("Indicate the threshold (\'d\' for default): ")
  temp <- as.character(readLines(n = 1))
  if(toupper(temp)!="D") { Dlim <- as.numeric(temp) }
}

## start the series if computations
cat(paste("Boundaries computed (required: ",ntry,")\n",sep=""))

## for loop
bdr.values <- -1 # used so that the first boundary is automatically the best

for(i in 0:(ntry-1)){
  temp <- monmonier(xy, dist, cn.nb,skip.local.diff=i,scanthres=FALSE,
                            threshold=Dlim, bd.length=bd.length, allowLoop=allowLoop)
  current.bdr.value <- sum(c(temp$run1$dir1$values, temp$run1$dir2$values))

  if(all(current.bdr.value > bdr.values)) {
      bdr.best <- temp
      bdr.skip <- i
  }

  bdr.values <- c(bdr.values , current.bdr.value)
  cat(paste(1+i," ")) # print progression in live
}

## remove the first value (-1) in bdr.value
bdr.values <- bdr.values[-1]

## graphical display
if(display.graph) barplot(bdr.values,xlab="Local differences skipped",ylab="Sum of all local differences",names.arg=0:(ntry-1))

## return the best value of skip.local.diff, or the corresponding object
if(!return.best) {
  return(bdr.skip)
  } else {
    cat(paste("\nOptimal number of skipped local differences: ",bdr.skip,"\n"))
    prevcall <- as.list(match.call())
    newcall <- bquote( monmonier(xy=.(prevcall$xy), dist=.(prevcall$dist), cn=.(prevcall$cn),
                                 threshold=.(bdr.best$threshold), bd.length=.(bd.length),
                                 nrun=1, skip.local.diff=.(as.numeric(bdr.skip)), scanthres=FALSE, allowLoop=.(allowLoop)) )

    ## assign the appropriate call to the result
    bdr.best$call <- newcall

    ## exp <- bquote( monmonier(xy=.(prevcall$xy),dist=.(prevcall$dist),cn=.(prevcall$cn),skip=.(bdr.skip),
    ##                          ,scan=FALSE,thres=.(Dlim),bd.length=.(bd.length),allowLoop=.(allowLoop)) )
    return(bdr.best)
    }
}
