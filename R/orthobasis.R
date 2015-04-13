## THIS CODE HAS BEEN TAKEN FROM ADE4
## IT HAS BEEN MOVED TO ADESPATIAL
## THIS IS MERELY A TEMPORARY FIX.

## only orthobasis.listw is used.
## renamed into .orthobasis.listw to avoid having to document it.

## - Thibaut Jombart, April 2015

#######################################################

## print.orthobasis <- function(x,...) {
##     if (!inherits(x,"orthobasis")) stop ("for 'orthobasis' object")
##     cat("Orthonormal basis: ")
##     n <- nrow(x)
##     p <- ncol(x)
##     if (n!=(p+1)) stop ("Non convenient dimension: author's error")
##     cat("data.frame with",n,"rows and",ncol(x),"columns\n")
##     cat("--------------------------------------\n")
##     cat("Columns are an orthonormal basis of 1n-orthogonal for\n")
##     cat("the inner product defined by the weights attribute\n")
##     cat("---------------------------------------\n")
##     w <- attributes(x)
##     if (!is.null(w$"names")) cat("names =", w$names[1],"...",w$names[p],"\n")
##     if (!is.null(w$"row.names")) cat("row.names =", w$row.names[1],"...",w$row.names[n],"\n")
##     if (!is.null(w$"weights")) cat("weights =", w$weights[1],"...",w$weights[n],"\n")
##     if (!is.null(w$"values")) cat("values =", w$values[1],"...",w$values[p],"\n")
##     if (!is.null(w$"class")) cat("class =", w$class,"\n")
##     if (!is.null(w$"call")) {
##         cat("call =")
##         print(w$"call")
##     }
##   }


## orthobasis.mat <- function(mat, cnw=TRUE) {
##     if (!is.matrix(mat)) stop ("matrix expected")
##     if (any(mat<0)) stop ("negative value in 'mat'")
##     if (nrow(mat)!=ncol(mat)) stop ("squared matrix expected")
##     mat <- (mat+t(mat))/2
##     nlig <- nrow(mat)
##     if (is.null(dimnames(mat))) {
##         w <- paste("P",1:nrow(mat),sep="")
##         dimnames(mat) <- list(w,w)
##     }
##     labels <- dimnames(mat)[[1]]
##     if (cnw) {
##         margi <- apply(mat,1,sum)
##         margi <- max(margi)-margi
##         mat <- mat+diag(margi)
##     }
##     mat <- mat/sum(mat)
##     wt <- rep ((1/nlig),nlig)
##     # calculs extensibles à une pondération quelconque
##     wt <- wt/sum(wt)
##     # si mat wt est la pondération marginale associée à mat
##     # tot = sum(mat)
##     # mat = mat-matrix(wt,nlig,nlig,byrow=TRUE)*wt*tot
##     # encore plus particulier mat = mat-1/nlig/nlig
##     # en général les précédents sont des cas particuliers
##     U <- matrix(1,nlig,nlig)
##     U  <- diag(1,nlig)-U*wt
##     mat <- U%*%mat%*%t(U)
##     wt <- sqrt(wt)
##     mat <- t(t(mat)/wt)
##     mat <- mat/wt
##     eig <- eigen(mat,sym=TRUE)
##     w0 <- abs(eig$values)/max(abs(eig$values))
##     tol <- 1e-07
##     w0 <- which(w0<tol)
##     if (length(w0)==0) stop ("abnormal output : no null eigenvalue")
##     else if (length(w0)==1) w0 <- (1:nlig)[-w0]
##     else if (length(w0)>1) {
##         # on ajoute le vecteur dérivé de 1n
##         w <- cbind(wt,eig$vectors[,w0])
##         # on orthonormalise l'ensemble
##         w <- qr.Q(qr(w))
##         # on met les valeurs propres à 0
##         eig$values[w0] <- 0
##         # on remplace les vecteurs du noyau par une base orthonormée contenant
##         # en première position le parasite
##         eig$vectors[,w0] <- w[,-ncol(w)]
##         # on enlève la position du parasite
##         w0 <- (1:nlig)[-w0[1]]
##     }
##     mat <- eig$vectors[,w0]/wt
##     mat <- data.frame(mat)
##     row.names(mat) <- labels
##     names(mat) <- paste("S",1:(nlig-1),sep="")
##     attr(mat,"values") <- eig$values[w0]
##     attr(mat,"weights") <- rep(1/nlig,nlig)
##     attr(mat,"call") <- match.call()
##     attr(mat,"class") <- c("orthobasis","data.frame")
##     return(mat)
## }

## "orthobasis.haar" <- function(n) {
## # on définit deux fonctions :
##     appel = match.call()
##     a <- log(n)/log(2)
##     b <- floor(a)
##     if ((a-b)^2>1e-10) stop ("Haar is not a power of 2")
## # la première est écrite par Daniel et elle donne la démonstration (par analogie avec la fonction qui construit la base Bscores)
## # que la base Bscores est exactement la base de Haar quand on prend une phylogénie régulière résolue.
## "haar.basis.1" <- function (n) {
##     pari <- matrix(c(1,n),1)
##     "div2" <- function (mat) {
##         res <- NULL
##         for (k in 1 : nrow(mat)) {
##             n1 <- mat[k,1]
##             n2 <- mat[k,2]
##             diff <- n2-n1
##             if (diff <=0) break
##             n3 <- floor((n1+n2)/2)
##             res <- rbind(res,c(n1,n3),c(n3+1,n2))
##         }
##         if (!is.null(res)) pari <<- rbind(pari,res)
##         return(res)
##     }
##     mat <- div2(pari)
##     while (!is.null(mat)) mat <- div2(mat)
##     res <- NULL
##     for (k in 1:nrow(pari)) {
##         x<-rep(0,n)
##         x[(pari[k,1]):(pari[k,2])] <- 1
##         res <-c(res,x)
##     }
##     res = matrix(res,n)
##     res <- qr.Q(qr(res))
##     res <- res[, -1] * sqrt(n)
##     res <- data.frame(res)
##     row.names(res) <- paste("u",1:n,sep="")
##     names(res) <- paste("B",1:(n-1),sep="")
## return(res)
## }

## # la seconde exploite les potentialités de la librairie waveslim, en remarquant qu'il existe un lien étroit entre la définition des filtres et la définition
## # des bases. Cette stratégie permettra à l'avenir de définir les bases associées à d'autres famille de fonctions.
## "haar.basis.2" <-  function (n) {
##     if (!require(waveslim)) stop ("Please install waveslim")
##     J <- a    #nombre de niveau
##     res <- matrix(0, nrow = n,ncol = n-1)
##     filter.seq <- "H" #filtre correspondant au niveau 1
##     h <- waveslim::wavelet.filter(wf.name = "haar", filter.seq = filter.seq)   #paramètre du filtre au niveau 1
##     k <- 0
##         for(i in 1:J){
##         z <- rep(h,2**(J-i))
##         x <- 1:n
##         y <- rep((n-1-k):(n-2**(J-i)-k),rep(2**i,2**(J-i)))
##         for(j in 1:n)   res[x[j],y[j]] <- z[j]
##         k <- k+2**(J-i)
##         filter.seq <- paste(filter.seq, "L", sep = "")
##         h <- waveslim::wavelet.filter(wf.name = "haar", filter.seq = filter.seq)
##         }
##         res <- res*sqrt(n)
##         res <- data.frame(res)
##         row.names(res) <- paste("u", 1:n, sep = "")
##         names(res) <- paste("B", 1:(n-1), sep = "")
## return(res)
## }

## # suivant que n est grand (n > 257) ou non, on choisit l'une des deux stratégies :
##     if (n < 257)
##         res <- haar.basis.1(n)
##         else
##             res <- haar.basis.2(n)

##     attr(res,"values") <- NULL
##     attr(res,"weights") <- rep(1/n,n)
##     attr(res,"call") <- appel
##     attr(res,"class") <- c("orthobasis","data.frame")
##     return(res)
## }

## "orthobasis.line" <- function (n) {
##     appel <- match.call()
##     # solution de Cornillon p. 12
##     res <- NULL
##     for (k in 1:(n-1)) {
##         x <- cos(k*pi*(2*(1:n)-1)/2/n)
##         x <- sqrt(n)*x/sqrt(sum(x*x))
##         res <-c(res,x)
##     }
##     res=matrix(res,n)
##     res <- data.frame(res)
##     row.names(res) <- paste("u",1:n,sep="")
##     names(res) <- paste("B",1:(n-1),sep="")
##     w <- (1:(n-1))*pi/2/n
##     valpro <- 4*(sin(w)^2)/n
##     poivoisi <- c(1,rep(2,n-2),1)
##     poivoisi <- poivoisi/sum(poivoisi)
##     norm <- unlist(apply(res, 2, function(a) sum(a*a*poivoisi)))
##     y <- valpro*n*n/2/(n-1)
##     val <- norm - y
##     attr(res,"values") <- val
##     attr(res,"weights") <- rep(1/n,n)
##     attr(res,"call") <- appel
##     attr(res,"class") <- c("orthobasis","data.frame")

##     # vérification locale. Ce paragraphe vérifie que les vecteurs et les valeurs
##     # proposée par Cornillon p. 12 sont bien les vecteurs propres de l'opérateur de voisinage
##     # rangée dans la solution analytique par variance locale croissante
##     # l'article de Méot est erroné et a donné le graphe circulaire pour le graphe linéaire
##     # d0=neig2mat(neig(n.lin=n))
##     # d0 = d0/n
##     # d1=apply(d0,1,sum)
##     # d0=diag(d1)-d0
##     # fun2 <- function(x) {
##     #     z <- sum(t(d0*x)*x)/n
##     #     z <- z/sum(x*x)
##     #     return(z)
##     # }
##     # lambda <- unlist(apply(res,2,fun2))
##     # print(lambda)
##     # print(attr(res,"values"))
##     # plot(lambda,attr(res,"values"))
##     # abline(lm(attr(res,"values")~lambda))
##     # print(coefficients(lm(attr(res,"values")~lambda)))

##     # vérification que les valeurs dérivées des valeurs propres sont exactement des indices de Moran
##     # d = neig2mat(neig(n.lin=n))
##     # d = d/sum(d) # Moran type W
##     # moran <- unlist(lapply(res,function(x) sum(t(d*x)*x)))
##     # print(moran)
##     # plot(moran,attr(res,"values"))
##     # abline(lm(attr(res,"values")~moran))
##     # print(summary(lm(attr(res,"values")~moran)))
##     return(res)
## }

## "orthobasis.circ" <- function (n) {
##     appel = match.call()
##     if (n<3) stop ("'n' too small")
##     "vecprosin" <- function(k) {
##         x <- sin(2*k*pi*(1:n)/n)
##         x <- x/sqrt(sum(x*x))
##     }
##     "vecprocos" <- function(k) {
##         x <- cos(2*k*pi*(1:n)/n)
##         x <- x/sqrt(sum(x*x))
##     }
##     "valpro" <- function(k,bis=TRUE) {
##         x <- (4/n)*((sin(k*pi/n))^2)
##         if (bis) x <- c(x,x)
##         return(x)
##     }

##     k <- floor(n/2)
##     if (k==n/2) {
##         #n est pair
##         w1 <- matrix(unlist(lapply(1:k,vecprocos)),n,k)
##         w2 <- matrix(unlist(lapply(1:(k-1),vecprosin)),n,k-1)
##         res <- cbind(w1,w2)
##         res[,seq(1,2*k-1,by=2)]<-w1
##         res[,seq(2,2*k-2,by=2)]<-w2
##         vp <- unlist(lapply(1:(k-1),valpro))
##         vp <- c(vp, valpro(k,FALSE))
##     } else {
##         # n est impair
##         w1 <- matrix(unlist(lapply(1:k,vecprocos)),n,k)
##         w2 <- matrix(unlist(lapply(1:k,vecprosin)),n,k)
##         res <- cbind(w1,w2)
##         res[,seq(1,2*k-1,by=2)]<-w1
##         res[,seq(2,2*k,by=2)]<-w2
##         vp <- unlist(lapply(1:k,valpro))
##     }
##     res=sqrt(n)*res
##     res <- as.data.frame(res)
##     row.names(res) <- paste("u",1:n,sep="")
##     names(res) <- paste("B",1:(n-1),sep="")
##     attr(res,"values") <- 1 - n*vp/2
##     attr(res,"weights") <- rep(1/n,n)
##     attr(res,"call") <- appel
##     attr(res,"class") <- c("orthobasis","data.frame")
##     # vérification qu'on a exactement des indices de Moran à partie des valeurs propres
##     # d = neig2mat(neig(n.cir=n))
##     # d = d/sum(d) # Moran type W
##     # moran <- unlist(lapply(res,function(x) sum(t(d*x)*x)))
##     # print(moran)
##     # plot(moran,attr(res,"values"))
##     # abline(lm(attr(res,"values")~moran))
##     # print(summary(lm(attr(res,"values")~moran)))
##     return(res)
## }

".orthobasis.listw" <- function( listw) {
    appel = match.call()
    if(!inherits(listw,"listw")) stop ("object of class 'listw' expected")
    if(listw$style!="W") stop ("object of class 'listw' with style 'W' expected")
    n = length(listw$weights)
    fun <- function (x) {
        num = listw$neighbours[[x]]
        wei = listw$weights[[x]]
        res = rep(0,n)
        res[num] = wei
        return (res)
    }
    b0 <- matrix(unlist(lapply(1:n,fun)),n,n)
    b0=(t(b0)+b0)/2
    b0=bicenter.wt(b0)
    a0 <- eigen(b0, sym = TRUE)
    #barplot(a0$values)
    a0 <- a0$vectors
    a0 <- cbind(rep(1,n),a0)
    a0 <- qr.Q(qr(a0))
    a0 <- as.data.frame(a0[,-1])*sqrt(n)
    row.names(a0) <- attr(listw,"region.id")
    names(a0) <- paste("VP", 1:(n-1), sep = "")
    z <- apply(a0,2,function(x) sum((t(b0*x)*x))/n)
    attr(a0,"values") <- z
    attr(a0,"weights") <- rep(1/n,n)
    attr(a0,"call") <- appel
    attr(a0,"class") <- c("orthobasis","data.frame")
    return(a0)
}


## "orthobasis.neig" <- function( neig) {
##     appel = match.call()
##     if(!inherits(neig,"neig")) stop ("object of class 'neig' expected")
##     n <- length(attr(neig,"degree"))
##     m <- sum(attr(neig,"degree"))
##     poivoisi <- attr(neig,"degree")/m
##     if (is.null(names(poivoisi))) names(poivoisi) <- as.character(1:n)
##     d0 = neig2mat(neig)
##     d0 = diag(poivoisi)-d0/m
##     eig <- eigen(d0, sym = TRUE)
##     ########
##     tol <- 1e-07
##     w0 <- abs(eig$values)/max(abs(eig$values))
##     w0 <- which(w0<tol)
##     if (length(w0)==0) stop ("abnormal output : no null eigenvalue")
##     else if (length(w0)==1) w0 <- (1:n)[-w0]
##     else if (length(w0)>1) {
##         # on ajoute le vecteur dérivé de 1n
##         wt <- rep(1,n)
##         w <- cbind(wt,eig$vectors[,w0])
##         # on orthonormalise l'ensemble
##         w <- qr.Q(qr(w))
##         # on met les valeurs propres à 0
##         eig$values[w0] <- 0
##         # on remplace les vecteurs du noyau par une base orthonormée contenant
##         # en première position le parasite
##         eig$vectors[,w0] <- w[,-ncol(w)]
##         # on enlève la position du parasite
##         w0 <- (1:n)[-w0[1]]
##     }
##     w0 <- rev(w0)
##     valpro <- eig$values[w0]
##     eig <- eig$vectors[,w0]
##     eig <- as.data.frame(eig)*sqrt(n)
##     z <- apply(eig,2,function(x) sum(x*x*poivoisi))
##     z <- z - valpro*n
##     w <- rev(order(z))
##     z <- z[w]
##     eig <- eig[,w]
##     row.names(eig) <- names(poivoisi)
##     names(eig) <- paste("VP", 1:(n-1), sep = "")
##     attr(eig,"values") <- z
##     attr(eig,"weights") <- rep(1/n,n)
##     attr(eig,"call") <- appel
##     attr(eig,"class") <- c("orthobasis","data.frame")
##     return(eig)
## }
