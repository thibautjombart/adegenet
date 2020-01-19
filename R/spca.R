###############################################
##
## spatial Principal Components Analysis
##
## require ade4, spdep and eventually tripack
##
## generic functions were derived from
## those of multispati class (ade4)
##
## T. Jombart (t.jombart@imperial.ac.uk)
## 31 may 2007
##############################################

spca <- function (...) UseMethod("spca")



################
## spca.default
################
#' @export
spca.default <- function(x, ...) {
  stop(sprintf("No spca method for object of class %s",
               paste(class(x), collapse = " ")))
}





###############
## spca.matrix
###############
#' @export
spca.matrix <- function(x, xy = NULL, cn = NULL, matWeight = NULL,
                        center = TRUE, scale = FALSE, scannf = TRUE,
                        nfposi = 1, nfnega = 1,
                        type = NULL, ask = TRUE,
                        plot.nb = TRUE, edit.nb = FALSE,
                        truenames = TRUE,
                        d1 = NULL, d2 = NULL, k = NULL,
                        a = NULL, dmin = NULL, ...) {

  ## check type of x: only numeric values are acceptable

  if (!is.numeric(x)) {
    stop("Only matrices of numeric values are accepted.")
  }


  ## check and handle xy coordinates

  if(is.null(xy) & (inherits(cn,"nb") & !inherits(cn,"listw")) ){
    xy <- attr(cn,"xy")  # xy can be retrieved from a nb object (not from listw)
  }
  if(is.null(xy)) {
    stop("xy coordinates are not provided")
  }
  if(is.data.frame(xy)) {
    xy <- as.matrix(xy)
  }
  if(!is.matrix(xy)) {
    stop("provided 'xy' cannot be converted to matrix")
  }
  if(ncol(xy) != 2) {
    stop("xy does not have two columns.")
  }
  if(nrow(xy) != nrow(x)) {
    stop("x and xy must have the same row numbers.")
  }


  ## find the spatial weights

  resCN <- NULL

  ## connection network from matWeight

  if (!is.null(matWeight)) {
    if (!is.matrix(matWeight)) {
      stop("matWeight is not a matrix")
    }
    if (!is.numeric(matWeight)) {
      stop("matWeight is not numeric")
    }
    if (nrow(matWeight) != ncol(matWeight)) {
      stop("matWeight is not square")
    }
    if (nrow(matWeight) != nrow(x)) {
      stop("dimension of datWeight does not match data")
    }
    diag(matWeight) <- 0
    matWeight <- prop.table(matWeight, 1)
    resCN <- mat2listw(matWeight)
    resCN$style <- "W"

  }

  ## connection network from cn argument

  if(is.null(resCN) & !is.null(cn)) {
    if(inherits(cn,"nb")) {
      if(!inherits(cn,"listw")){ # cn is a 'pure' nb object (i.e., nb but not listw)
        cn <- nb2listw(cn, style="W", zero.policy=TRUE)
      }
      resCN <- cn
    } else {
      stop("cn does not have a recognized class")
    }
  }


  ## connection network from xy coordinates
  if(is.null(resCN)) {
    resCN <- chooseCN(xy=xy, ask=ask, type=type, plot.nb=plot.nb, edit.nb=edit.nb,
                      result.type="listw", d1=d1, d2=d2, k=k, a=a, dmin=dmin)
  }


  ## perform the analyses: basic PCA followed by multispati

  x_pca <- ade4::dudi.pca(x, center = center, scale = scale, scannf = FALSE)

  out <- ade4::multispati(dudi = x_pca, listw = resCN, scannf = scannf,
                          nfposi = nfposi, nfnega = nfnega)

  nfposi <- out$nfposi
  nfnega <- out$nfnega

  out$tab <- x_pca$tab
  out$xy <- xy
  rownames(out$xy) <- rownames(out$li)
  colnames(out$xy) <- c("x","y")

  out$lw <- resCN

  dots <- list(...)
  if (!is.null(dots$call)) {
    out$call <- dots$call
  } else {
    out$call <- match.call()
  }

  posaxes <- if (nfposi > 0) {1:nfposi} else NULL
  negaxes <- if (nfnega > 0) {(length(out$eig)-nfnega+1):length(out$eig)} else NULL
  keptaxes <- c(posaxes, negaxes)

  ## set names of different components
  colnames(out$c1) <- paste("Axis",keptaxes)
  colnames(out$li) <- paste("Axis",keptaxes)
  colnames(out$ls) <- paste("Axis",keptaxes)
  row.names(out$c1) <- colnames(x)
  colnames(out$as) <- colnames(out$c1)
  temp <- row.names(out$as)
  row.names(out$as) <- paste("PCA", temp)

  class(out) <- "spca"

  return(out)

}






###################
## spca.data.frame
###################

#' @method spca data.frame
#' @export
spca.data.frame <- function(x, xy = NULL, cn = NULL, matWeight = NULL,
                            center = TRUE, scale = FALSE, scannf = TRUE,
                            nfposi = 1, nfnega = 1,
                            type = NULL, ask = TRUE,
                            plot.nb = TRUE, edit.nb = FALSE,
                            truenames = TRUE,
                            d1 = NULL, d2 = NULL, k = NULL,
                            a = NULL, dmin = NULL, ...) {

  call <- match.call()

  spca(as.matrix(x), xy = xy, cn = cn, matWeight = matWeight, center = center,
       cale = scale, scannf = scannf, nfposi = nfposi, nfnega = nfnega,
       type = type, ask = ask, plot.nb = plot.nb, edit.nb = edit.nb,
       truenames = truenames, d1 = d1, d2 = d2, k = k, a = a, dmin = dmin,
       call = call, ...)
}






################
## spca genind
################

#' @export
spca.genind <- function(obj, xy = NULL, cn = NULL, matWeight = NULL,
                        scale = FALSE, scannf = TRUE,
                        nfposi = 1, nfnega = 1,
                        type = NULL, ask = TRUE,
                        plot.nb = TRUE, edit.nb = FALSE,
                        truenames = TRUE,
                        d1 = NULL, d2 = NULL, k = NULL,
                        a = NULL, dmin = NULL, ...){

  ## first checks

  invisible(validObject(obj))


  ## handle xy coordinates
  if(is.null(xy) & !is.null(obj$other$xy)) {
    xy <- obj$other$xy # xy from @other$xy if it exists
  }



  ## == spatial weights are done ==


  ## handle NAs warning
  if(any(is.na(obj@tab))){
    warning("NAs in data are automatically replaced (to mean allele frequency)")
  }

  ## handle NAs, centring and scaling
  X <- tab(obj, freq = TRUE, NA.method = "mean")

  call <- match.call()

  spca(X, xy = xy, cn = cn, matWeight = matWeight,
       center = TRUE, scale = scale, scannf = scannf,
       nfposi = nfposi, nfnega = nfnega,
       type = type, ask = ask,
       plot.nb = plot.nb, edit.nb = edit.nb,
       truenames = truenames,
       d1 = d1, d2 = d2, k = k,
       a = a, dmin = dmin,
       call = call, ...)

} # end spca.genind






################
## spca genpop
################

#' @export
spca.genpop <- function(obj, xy = NULL, cn = NULL, matWeight = NULL,
                        scale = FALSE, scannf = TRUE,
                        nfposi = 1, nfnega = 1,
                        type = NULL, ask = TRUE,
                        plot.nb = TRUE, edit.nb = FALSE,
                        truenames = TRUE,
                        d1 = NULL, d2 = NULL, k = NULL,
                        a = NULL, dmin = NULL, ...){

  ## first checks

  invisible(validObject(obj))


  ## handle xy coordinates
  if(is.null(xy) & !is.null(obj$other$xy)) {
    xy <- obj$other$xy # xy from @other$xy if it exists
  }


  ## handle NAs warning
  if(any(is.na(obj@tab))){
    warning("NAs in data are automatically replaced (to mean allele frequency)")
  }

  ## handle NAs, centring and scaling
  X <- tab(obj, freq = TRUE, NA.method = "mean")

  call <- match.call()

  spca(X, xy = xy, cn = cn, matWeight = matWeight,
       center = TRUE, scale = scale, scannf = scannf,
       nfposi = nfposi, nfnega = nfnega,
       type = type, ask = ask,
       plot.nb = plot.nb, edit.nb = edit.nb,
       truenames = truenames,
       d1 = d1, d2 = d2, k = k,
       a = a, dmin = dmin,
       call = call, ...)

} # end spca.genpop






######################
## Function print.spca
######################

#' @method print spca
#' @export
print.spca <- function(x, ...){
  cat("\t########################################\n")
  cat("\t# spatial Principal Component Analysis #\n")
  cat("\t########################################\n")
  cat("class: ")
  cat(class(x))
  cat("\n$call: ")
  print(x$call)
  cat("\n$nfposi:", x$nfposi, "axis-components saved")
  cat("\n$nfnega:", x$nfnega, "axis-components saved")

  cat("\nPositive eigenvalues: ")
  l0 <- sum(x$eig >= 0)
  cat(signif(x$eig, 4)[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...\n")
  else cat("\n")
  cat("Negative eigenvalues: ")
  l0 <- sum(x$eig <= 0)
  cat(sort(signif(x$eig, 4))[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...\n")
  else cat("\n")
  cat('\n')
  sumry <- array("", c(1, 4), list(1, c("vector", "length",
                                        "mode", "content")))
  sumry[1, ] <- c('$eig', length(x$eig), mode(x$eig), 'eigenvalues')
  class(sumry) <- "table"
  print(sumry)
  cat("\n")
  sumry <- array("", c(5, 4),
                 list(1:5, c("data.frame", "nrow", "ncol", "content")))
  sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab),
                  "transformed data: optionally centred / scaled")
  sumry[2, ] <- c("$c1", nrow(x$c1), ncol(x$c1),
                  "principal axes: scaled vectors of alleles loadings")
  sumry[3, ] <- c("$li", nrow(x$li), ncol(x$li),
                  "principal components: coordinates of entities ('scores')")
  sumry[4, ] <- c("$ls", nrow(x$ls), ncol(x$ls),
                  "lag vector of principal components")
  sumry[5, ] <- c("$as", nrow(x$as), ncol(x$as),
                  "pca axes onto spca axes")

  class(sumry) <- "table"
  print(sumry)

  cat("\n$xy: matrix of spatial coordinates")
  cat("\n$lw: a list of spatial weights (class 'listw')")

  cat("\n\nother elements: ")
  if (length(names(x)) > 10)
    cat(names(x)[11:(length(names(x)))], "\n")
  else cat("NULL\n")
}





########################
## Function summary.spca
########################
#' @method summary spca
#' @export
summary.spca <- function (object, ..., printres=TRUE) {
  if (!inherits(object, "spca"))stop("to be used with 'spca' object")

  #util <- function(n) { ## no longer used
  #  x <- "1"
  #  for (i in 2:n) x[i] <- paste(x[i - 1], i, sep = "+")
  #  return(x)
  #}
  norm.w <- function(X, w) {
    f2 <- function(v) sum(v * v * w)/sum(w)
    norm <- apply(X, 2, f2)
    return(norm)
  }

  resfin <- list()

  if(printres) {
    cat("\nSpatial principal component analysis\n")
    cat("\nCall: ")
    print(object$call)
  }

  appel <- as.list(object$call)
  ## compute original pca
  # prepare data
  obj <- eval(appel$obj)
  if(is.null(appel$truenames)) appel$truenames <- FALSE

  f1 <- function(vec){
    m <- mean(vec,na.rm=TRUE)
    vec[is.na(vec)] <- m
    return(vec)
  }

  if(is.genind(obj)) { X <- obj@tab }
  if(is.genpop(obj)) { X <- makefreq(obj, quiet=TRUE) }

  X <- apply(X,2,f1)

  nfposi <- object$nfposi
  nfnega <- object$nfnega

  dudi <- dudi.pca(X, center=TRUE, scale=FALSE, scannf=FALSE, nf=nfposi+nfnega)
  ## end of pca

  lw <- object$lw

  # I0, Imin, Imax
  n <- nrow(X)
  I0 <- -1/(n-1)
  L <- listw2mat(lw)
  ## use 'as.numeric' to avoid possible bug with large matrices,
  ## returning complex numbers with zero imaginary parts
  eigL <- suppressWarnings(as.numeric(eigen(0.5*(L+t(L)))$values))
  Imin <- min(eigL)
  Imax <- max(eigL)
  Ival <- data.frame(I0=I0,Imin=Imin,Imax=Imax)
  row.names(Ival) <- ""
  if(printres) {
    cat("\nConnection network statistics:\n")
    print(Ival)
  }

  Istat <- c(I0,Imin,Imax)
  names(Istat) <- c("I0","Imin","Imax")
  resfin$Istat <- Istat


  # les scores de l'analyse de base
  nf <- dudi$nf
  eig <- dudi$eig[1:nf]
  cum <- cumsum(dudi$eig)[1:nf]
  ratio <- cum/sum(dudi$eig)
  w <- apply(dudi$l1,2,lag.listw,x=lw)
  moran <- apply(w*as.matrix(dudi$l1)*dudi$lw,2,sum)
  res <- data.frame(var=eig,cum=cum,ratio=ratio, moran=moran)
  row.names(res) <- paste("Axis",1:nf)
  if(printres) {
    cat("\nScores from the centred PCA\n")
    print(res)
  }

  resfin$pca <- res


  # les scores de l'analyse spatiale
  # on recalcule l'objet en gardant tous les axes
  eig <- object$eig
  nfposimax <- sum(eig > 0)
  nfnegamax <- sum(eig < 0)

  ms <- multispati(dudi=dudi, listw=lw, scannf=FALSE,
                   nfposi=nfposimax, nfnega=nfnegamax)

  ndim <- dudi$rank
  nf <- nfposi + nfnega
  agarder <- c(1:nfposi,if (nfnega>0) (ndim-nfnega+1):ndim)
  varspa <- norm.w(ms$li,dudi$lw)
  moran <- apply(as.matrix(ms$li)*as.matrix(ms$ls)*dudi$lw,2,sum)
  res <- data.frame(eig=eig,var=varspa,moran=moran/varspa)
  row.names(res) <- paste("Axis",1:length(eig))

  if(printres) {
    cat("\nsPCA eigenvalues decomposition:\n")
    print(res[agarder,])
  }

  resfin$spca <- res

  return(invisible(resfin))
}



#####################
## Function plot.spca
#####################
#' @method plot spca
#' @export
plot.spca <- function (x, axis = 1, useLag=FALSE, ...){
  if (!inherits(x, "spca")) stop("Use only with 'spca' objects.")

  if(axis>ncol(x$li)) stop("wrong axis required.")

  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mar = rep(.1,4), mfrow=c(3,2))

  n <- nrow(x$li)
  xy <- x$xy

  ## handle useLag argument
  if(useLag){
    z <- x$ls[,axis]
  } else {
    z <- x$li[,axis]
  } # end if useLag
  nfposi <- x$nfposi
  nfnega <- x$nfnega
  ## handle neig parameter - hide cn if nore than 100 links
  nLinks <- sum(card(x$lw$neighbours))
  if(nLinks < 500) {
    neig <- nb2neig(x$lw$neighbours)
  } else {
    neig <- NULL
  }

  sub <- paste("Score",axis)
  csub <- 2

  # 1
  if(n<30) clab <- 1 else clab <- 0
  s.label(xy, clabel=clab, include.origin=FALSE, addaxes=FALSE, neig=neig,
          cneig=1, sub="Connection network", csub=2)

  # 2
  s.image(xy,z, include.origin=FALSE, grid=TRUE, kgrid=10, cgrid=1,
          sub=sub, csub=csub, possub="bottomleft")
  box()

  # 3
  if(n<30) {neig <- nb2neig(x$lw$neighbours)} else {neig <- NULL}
  s.value(xy,z, include.origin=FALSE, addaxes=FALSE, clegend=0, csize=.6,
          neig=neig, sub=sub, csub=csub, possub="bottomleft")

  # 4
  s.value(xy,z, include.origin=FALSE, addaxes=FALSE, clegend=0, csize=.6,
          method="greylevel", neig=neig, sub=sub, csub=csub, possub="bottomleft")

  # 5
  omar <- par("mar")
  par(mar = c(0.8, 2.8, 0.8, 0.8))
  m <- length(x$eig)
  col.w <- rep("white", m) # elles sont toutes blanches
  col.w[1:nfposi] <- "grey"
  if (nfnega>0) {col.w[m:(m-nfnega+1)] <- "grey"}
  j <- axis
  if (j>nfposi) {j <- j-nfposi +m -nfnega}
  col.w[j] <- "black"
  barplot(x$eig, col = col.w)
  scatterutil.sub(cha ="Eigenvalues", csub = 2.5, possub = "topright")
  par(mar=rep(.1,4))
  box()
  par(mar=omar)

  # 6
  par(mar=c(4,4,2,1))
  screeplot(x,main="Eigenvalues decomposition")
  par(mar=rep(.1,4))
  box()
  return(invisible(match.call()))
}



##########################
## Function screeplot.spca
##########################
#' @method screeplot spca
#' @export
screeplot.spca <- function(x,...,main=NULL){

  opar <- par("las")
  on.exit(par(las=opar))

  sumry <- summary(x,printres=FALSE)

  labels <- lapply(1:length(x$eig),function(i) bquote(lambda[.(i)]))

  par(las=1)

  xmax <- sumry$pca[1,1]*1.1
  I0 <- sumry$Istat[1]
  Imin <- sumry$Istat[2]
  Imax <- sumry$Istat[3]

  plot(x=sumry$spca[,2],y=sumry$spca[,3],type='n',xlab='Variance',ylab="Spatial autocorrelation (I)",xlim=c(0,xmax),ylim=c(Imin*1.1,Imax*1.1),yaxt='n',...)
  text(x=sumry$spca[,2],y=sumry$spca[,3],do.call(expression,labels))

  ytick <- c(I0,round(seq(Imin,Imax,le=5),1))
  ytlab <- as.character(round(seq(Imin,Imax,le=5),1))
  ytlab <- c(as.character(round(I0,1)),as.character(round(Imin,1)),ytlab[2:4],as.character(round(Imax,1)))
  axis(side=2,at=ytick,labels=ytlab)

  rect(0,Imin,xmax,Imax,lty=2)
  segments(0,I0,xmax,I0,lty=2)
  abline(v=0)

  if(is.null(main)) main <- ("Spatial and variance components of the eigenvalues")
  title(main)

  return(invisible(match.call()))
}





###################
## colorplot method
###################
#' @method colorplot spca
#' @export
colorplot.spca <- function(x, axes=1:ncol(x$li), useLag=FALSE, ...){
  ## some checks
  if(!any(inherits(x,"spca"))) stop("x in not a spca object.")

  ## get args to be passed to colorplot
  xy <- x$xy

  if(useLag) {
    X <- as.matrix(x$ls)
  } else {
    X <- as.matrix(x$li)
  }

  ## call to colorplot
  colorplot(xy, X, axes, ...)

} # end colorplot.spca
