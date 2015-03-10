###################################
#
# Distances functions
# based on old ade4 'dist.genet'
# by Daniel Chessel
#
# Thibaut Jombart
# t.jombart@imperial.ac.uk
###################################



############################
# S3 method dist for genpop
############################
dist.genpop <- function(x, method = 1, diag = FALSE, upper = FALSE) {

  if(!is.genpop(x)) stop("x is not a valid genpop object")

  ## haploidy kludge (have to get rid of that later)
  if(x@ploidy==as.integer(1)){
  x@tab <- x@tab * 2
  x@ploidy <- as.integer(2)
  }


  ## check marker type
  checkType(x)


  METHODS = c("Nei","Edwards","Reynolds","Rodgers","Provesti")
  if (all((1:5)!=method)) {
    cat("1 = Nei 1972\n")
    cat("2 = Edwards 1971\n")
    cat("3 = Reynolds, Weir and Coockerman 1983\n")
    cat("4 = Rodgers 1972\n")
    cat("5 = Provesti 1975\n")
    cat("Select an integer (1-5): ")
    method <- as.integer(readLines(n = 1))
  }
  if (all((1:5)!=method)) (stop ("Non convenient method number"))

  nloc <- length(levels(x@loc.fac))
  loc.fac <- x@loc.fac
  X <- makefreq(x,missing="mean",quiet=TRUE)$tab
  # X is a matrix of allelic frequencies
  nlig <- nrow(X)

  if (method == 1) { # Nei
    d <- X%*%t(X)
    vec <- sqrt(diag(d))
    d <- d/vec[col(d)]
    d <- d/vec[row(d)]
    d <- -log(d)
    d <- as.dist(d)
  } else if (method == 2) { # Edward's (angular)
        X <- sqrt(X)
        d <- X%*%t(X)
        d <- 1-d/nloc
        diag(d) <- 0
        d <- sqrt(d)
        d <- as.dist(d)
    } else if (method == 3) { # Coancestrality coef (Reynold's)
       denomi <- X%*%t(X)
       vec <- apply(X,1,function(x) sum(x*x))
       d <- -2*denomi + vec[col(denomi)] + vec[row(denomi)]
       diag(d) <- 0
       denomi <- 2*nloc - 2*denomi
       diag(denomi) <- 1
       d <- d/denomi
       d <- sqrt(d)
       d <- as.dist(d)
    } else if (method == 4) { # Rogers' distance
      # kX is a list of K=nloc matrices
      kX <- lapply(split(X,loc.fac[col(X)]),matrix,nrow=nlig)
      dcano <- function(mat) {
        daux <- mat%*%t(mat)
        vec <- diag(daux)
        daux <- -2*daux + vec[col(daux)] + vec[row(daux)]
        diag(daux) <- 0
        daux <- sqrt(.5*daux)
        return(daux)
      }

      d <- matrix(0,nlig,nlig)
      for(i in 1:length(kX)) {
        d <- d + dcano(kX[[i]])
      }
      d <- d/length(kX)
      d <- as.dist(d)
    } else if (method ==5) { # Provesti (absolute genetic distance)
      w0 <- 1:(nlig-1)
      loca <- function(k) {
        w1 <- (k+1):nlig
        resloc <- unlist(lapply(w1, function(x) sum(abs(X[k,]-X[x,]))))
        return(resloc/(2*nloc))
      }
      d <- unlist(lapply(w0,loca))
    }
    attr(d, "Size") <- nlig
    attr(d, "Labels") <- x@pop.names
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- METHODS[method]
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)

} # end method dist for genpop
