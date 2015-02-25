#' Returns original points in results paths of an object of class 'monmonier'
#' 
#' The original implementation of \code{monmonier} in package \bold{adegenet}
#' returns path coordinates, \code{coords.monmonier} additionally displays
#' identities of the original points of the network, based on original
#' coordinates.
#' 
#' 
#' @param x an object of class \code{\link{monmonier}}.
#' @return Returns a list with elements according to the \code{x$nrun} result
#' of the \code{\link{monmonier}} object.  Corresponding path points are in the
#' same order as in the original object.
#' 
#' run1 (run2, ...): for each run, a list containing a matrix giving the
#' original points in the network (\code{first} and \code{second}, indicating
#' pairs of neighbours). Path coordinates are stored in columns \code{x.hw} and
#' \code{y.hw}. \code{first} and \code{second} are integers referring to the
#' row numbers in the \code{x$xy} matrix of the original
#' \code{\link{monmonier}} object.
#' @author Peter Solymos, \email{Solymos.Peter@@aotk.szie.hu},
#' \url{http://www.univet.hu/users/psolymos/personal/}
#' @seealso \code{\link{monmonier}}
#' @keywords methods manip
#' @examples
#' 
#' 
#' \dontrun{
#' if(require(spdep)){
#' 
#' load(system.file("files/mondata1.rda",package="adegenet"))
#' cn1 <- chooseCN(mondata1$xy,type=2,ask=FALSE)
#' mon1 <- monmonier(mondata1$xy,dist(mondata1$x1),cn1,threshold=2,nrun=3)
#' 
#' mon1$run1
#' mon1$run2
#' mon1$run3
#' path.coords <- coords.monmonier(mon1)
#' path.coords
#' }
#' }
#' 
#' @export coords.monmonier
coords.monmonier <- function(x){

# require(adegenet) no longer use since in adegenet
if (!inherits(x, "monmonier")) stop("Use only with 'monmonier' objects")

xy.full <- x$xy
n.points <- nrow(xy.full)
k <- 1

#cn list to cn matrix
cn.matr <- matrix(data = 0, nrow = n.points, ncol = n.points)

for(i in 1:n.points){
	eval <- is.element(c(1:n.points), x$cn[[i]])
	cn.matr[i,which(eval == TRUE)] <- 1
	}

#halfway <- matrix(data = NA, nrow = (n.points^2-n.points)/2, ncol = 4)
halfway <- matrix(data = NA, nrow = sum(cn.matr)/2, ncol = 4)
colnames(halfway) <- c("x.hw","y.hw","first","second")

for(first in 1:(n.points-1)){
for(second in (first+1):n.points){
	if(cn.matr[first,second] == 1){
	halfway[k,] <- c(
		(xy.full[first,1]+xy.full[second,1])/2,
		(xy.full[first,2]+xy.full[second,2])/2,
		first, second)
	k <- k+1}
	}}

output=list()
for(run in 1:x$nrun){
	runname <- paste('run',run,sep='')
	output[[runname]] <- list(dir1=list(),dir2=list())

dir1.in <- x[[runname]]$dir1$path
output[[runname]]$dir1 <- matrix(data = NA, nrow = nrow(dir1.in), ncol = 4)
colnames(output[[runname]]$dir1) <- c("x.hw","y.hw","first","second")
rownames(output[[runname]]$dir1) <- rownames(x[[runname]]$dir1$path)
for(i in 1:nrow(dir1.in)){
	eval.x <- is.element(halfway[,1], dir1.in[i,1])
	eval.y <- is.element(halfway[,2], dir1.in[i,2])
	output[[runname]]$dir1[i,] <- halfway[which(eval.x == TRUE & eval.y == TRUE),]
	}

dir2.in <- x[[runname]]$dir2$path
output[[runname]]$dir2 <- matrix(data = NA, nrow = nrow(dir2.in), ncol = 4)
colnames(output[[runname]]$dir2) <- c("x.hw","y.hw","first","second")
rownames(output[[runname]]$dir2) <- rownames(x[[runname]]$dir2$path)
for(i in 1:nrow(dir2.in)){
	eval.x <- is.element(halfway[,1], dir2.in[i,1])
	eval.y <- is.element(halfway[,2], dir2.in[i,2])
	output[[runname]]$dir2[i,] <- halfway[which(eval.x == TRUE & eval.y == TRUE),]
	}
}

return(output)
}
