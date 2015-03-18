#########################
# Function genind2genpop
#########################


#' Conversion from a genind to a genpop object
#'
#' The function \code{genind2genpop} converts genotypes data (genind) into
#' alleles counts per population (genpop).
#'
#' === 'missing' argument ===\cr The values of the 'missing' argument in
#' \code{genind2genpop} have the following effects:\cr - "NA": if all genotypes
#' of a population for a given allele are missing, count value will be NA\cr
#'
#' - "0": if all genotypes of a population for a given allele are missing,
#' count value will be 0\cr
#'
#' - "chi2": if all genotypes of a population for a given allele are missing,
#' count value will be that of a theoretical count in of a Chi-squared test.
#' This is obtained by the product of the margins sums divided by the total
#' number of alleles.\cr
#'
#' === processing the \code{@@other} slot ===\cr Essentially,
#' \code{genind2genpop} is about aggregating data per population. The function
#' can do the same for all numeric items in the \code{@@other} slot provided
#' they have the same length (for vectors) or the same number of rows
#' (matrix-like objects) as the number of genotypes. When the case is
#' encountered and if \code{process.other} is TRUE, then these objects are
#' processed using the function defined in \code{other.action} per population.
#' For instance, spatial coordinates of genotypes would be averaged to obtain
#' population coordinates.
#'
#' @param x an object of class \code{genind}.
#' @param pop a factor giving the population of each genotype in 'x'. If note
#' provided, sought in x@@pop, but if given, the argument prevails on x@@pop.
#' @param missing can be "NA", "0", or "chi2". See details for more
#' information.
#' @param quiet logical stating whether a conversion message must be printed
#' (TRUE,default) or not (FALSE).
#' @param process.other a logical indicating whether the \code{@@other} slot
#' should be processed (see details).
#' @param other.action a function to be used when processing the \code{@@other}
#' slot. By default, 'mean' is used.
#' @return A genpop object. The component @@other in 'x' is passed to the
#' created genpop object.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \linkS4class{genind}, \linkS4class{genpop},
#' \code{\link{na.replace}}
#' @keywords classes manip multivariate
#' @examples
#'
#' ## simple conversion
#' data(nancycats)
#' nancycats
#' catpop <- genind2genpop(nancycats)
#' catpop
#' summary(catpop)
#'
#' ## processing the @@other slot
#' data(sim2pop)
#' sim2pop$other$foo <- letters
#' sim2pop
#' dim(sim2pop$other$xy) # matches the number of genotypes
#' sim2pop$other$foo # does not match the number of genotypes
#'
#' obj <- genind2genpop(sim2pop, process.other=TRUE)
#' obj$other # the new xy is the populations' centre
#'
#' pch <- as.numeric(pop(sim2pop))
#' col <- pop(sim2pop)
#' levels(col) <- c("blue","red")
#' col <- as.character(col)
#' plot(sim2pop$other$xy, pch=pch, col=col)
#' text(obj$other$xy, lab=row.names(obj$other$xy), col=c("blue","red"), cex=2, font=2)
#'
#'
#' @export genind2genpop
genind2genpop <- function(x,pop=NULL,missing=c("NA","0","chi2"),quiet=FALSE,
                          process.other=FALSE, other.action=mean){

  if(!is.genind(x)) stop("x is not a valid genind object")
  checkType(x)

  if(is.null(x@pop) && is.null(pop)) {
      if(!quiet) warning("\npop is not provided either in x or in pop - assuming one single group")
      pop <- factor(rep(1, nrow(x@tab)))
  }

  missing <- match.arg(missing)

  if(!quiet) cat("\n Converting data from a genind to a genpop object... \n")

  ## choose pop argument over x@pop
   if(!is.null(pop)) {
    if(length(pop) != nrow(x@tab)) stop("inconsistent length for factor pop")
    # keep levels in order of appearance
    pop <- as.character(pop)
    pop <- factor(pop, levels=unique(pop))
  } else {
    pop <- pop(x)
    # keep levels in order of appearance
    ##pop <- as.character(pop)
    ##pop <- factor(pop, levels=unique(pop))
    ##if(!is.null(x@pop.names)) levels(pop) <- x@pop.names # restore real names
  }

  # make generic pop labels, store real pop names
  # pop.names <- levels(pop) ## no longer used

  # tabcount is a matrix pop x alleles, counting alleles per pop
  # *ploidy to have alleles counts
  f1 <- function(v){
    if(all(is.na(v))) return(NA) else return(sum(v,na.rm=TRUE))
  }

  f2 <- function(v){
    if(all(is.na(v)) || sum(v,na.rm=TRUE)==0) return(NA)
    return(v/(sum(v,na.rm=TRUE)))
  }

  tabcount <- apply(x@tab,2,function(c) tapply(c,pop,f1))
  tabcount <- round(tabcount,digits=0)
  # restitute matrix class when only one pop
  if(is.null(dim(tabcount))) {
    lab.col <- names(tabcount)
    tabcount <- matrix(tabcount,nrow=1)
    colnames(tabcount) <- lab.col
    rownames(tabcount) <- levels(pop)[1]
  }

  ## make final object
  if(x@type=="codom"){
      temp <- paste(rep(x@loc.names,x@loc.nall),unlist(x@all.names),sep=".")
  } else{
      temp <- x@loc.names
  }
  colnames(tabcount) <- temp

  prevcall <- match.call()

  res <- genpop(tab=tabcount, prevcall=prevcall, ploidy=x@ploidy, type=x@type)

  ## handle @other here
  res@other <- x@other
  if(process.other){
      ## auxiliary function doing the job
      fOther <- function(e){
          N <- nrow(x@tab)
          if(is.vector(e) && is.numeric(e) && length(e)==N){ # numeric vector
              res <- tapply(e, pop, other.action)
              return(res)
          } else if(is.matrix(e) && is.numeric(e) && nrow(e)==N){ # numeric matrix
              res <- apply(e, 2, function(vec) tapply(vec, pop, other.action))
              colnames(res) <- colnames(e)
              return(res)
          } else if(is.data.frame(e) && nrow(e)==N && all(sapply(e,is.numeric)) ){ # df of numeric vectors
              res <- lapply(e, function(vec) tapply(vec, pop, other.action))
              res <- data.frame(res)
              names(res) <- names(e)
              return(res)
          } else return(e)
      } # end fOther

      res@other <- lapply(res@other, fOther)
  } # end if(process.other)

  if(missing != "NA"){
      res <- na.replace(res, method=missing, quiet=quiet)
  }

  if(!quiet) cat("\n...done.\n\n")

  return(res)

} # end genind2genpop
