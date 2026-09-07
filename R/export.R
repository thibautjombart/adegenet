############################################
#
# Functions to transform a genind object
# into other R classes
#
# Thibaut Jombart
# t.jombart@imperial.ac.uk
#
############################################


## ############################
## # Function genind2hierfstat
## ############################
## genind2hierfstat <- function(x,pop=NULL){
##     ##  if(!inherits(x,"genind")) stop("x must be a genind object (see ?genind)")
##     ##   invisible(validObject(x))
##     if(!is.genind(x)) stop("x is not a valid genind object")
##     if(any(ploidy(x) != 2L)) stop("not implemented for non-diploid genotypes")
##     checkType(x)

##     if(is.null(pop)) pop <- pop(x)
##     if(is.null(pop)) pop <- as.factor(rep("P1",nrow(x@tab)))

##     ## ## NOTES ON THE CODING IN HIERFSTAT ##
##     ## - interpreting function is genot2al
##     ## - same coding has to be used for all loci
##     ## (i.e., all based on the maximum number of digits to be used)
##     ## - alleles have to be coded as integers
##     ## - alleles have to be sorted by increasing order when coding a genotype
##     ## - for instance, 121 is 1/21, 101 is 1/1, 11 is 1/1

##     ## find max number of alleles ##
##     max.nall <- max(nAll(x))
##     x@all.names <- lapply(alleles(x), function(e) .genlab("",max.nall)[1:length(e)])


##     ## VERSION USING GENIND2DF ##
##     gen <- genind2df(x, sep="", usepop=FALSE)
##     gen <- as.matrix(data.frame(lapply(gen, as.numeric)))
##     res <- cbind(as.numeric(pop),as.data.frame(gen))
##     colnames(res) <- c("pop", locNames(x))
##     if(!any(table(indNames(x))>1)){
##         rownames(res) <- indNames(x)
##     } else {
##         warning("non-unique labels for individuals; using generic labels")
##         rownames(res) <- 1:nrow(res)
##     }

##     return(res)
## } # end genind2hierfstat





#####################
# Function genind2df
#####################
#' Convert a genind object to a data.frame.
#'
#' The function \code{genind2df} converts a \linkS4class{genind} back to a
#' data.frame of raw allelic data.
#'
#' @aliases genind2df
#'
#' @param x a \linkS4class{genind} object
#' @param pop an optional factor giving the population of each individual.
#' @param sep a character string separating alleles. See details.
#' @param usepop a logical stating whether the population (argument \code{pop}
#' or \code{x@@pop} should be used (TRUE, default) or not (FALSE)).
#' @param oneColPerAll a logical stating whether or not alleles should be split
#' into columns (defaults to \code{FALSE}). This will only work with data with
#' consistent ploidies.
#'
#' @details With \code{oneColPerAll = TRUE}, observed allele copies are
#' retained and unobserved copies are padded with real \code{NA} values.
#' Zero-copy genotypes or genotypes containing missing allele counts produce
#' missing values in all allele columns. Allele order does not imply phase.
#' Genotypes with more copies than the declared ploidy are rejected.
#' This handling of incomplete genotypes follows the problem and approach
#' described by \code{ntakebay} in GitHub issue 320.
#'
#' @return a data.frame of raw allelic data, with individuals in rows and loci in column
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @seealso \code{\link{df2genind}}, \code{\link{import2genind}}, \code{\link{read.genetix}},
#' \code{\link{read.fstat}}, \code{\link{read.structure}}
#' @keywords manip
#' @examples
#'
#' ## simple example
#' df <- data.frame(locusA=c("11","11","12","32"),
#' locusB=c(NA,"34","55","15"),locusC=c("22","22","21","22"))
#' row.names(df) <- .genlab("genotype",4)
#' df
#'
#' obj <- df2genind(df, ploidy=2, ncode=1)
#' obj
#' obj@@tab
#'
#'
#' ## converting a genind as data.frame
#' genind2df(obj)
#' genind2df(obj, sep="/")
#'
#' @export
#'
genind2df <- function(x, pop=NULL, sep="", usepop=TRUE, oneColPerAll = FALSE){

  if(!is.genind(x)) stop("x is not a valid genind object")
  ## checkType(x)

  if(is.null(pop)) {
      pop <- x@pop
  }

  ## PA case ##
  if(x@type=="PA"){
      res <- tab(x)
      if(usepop && !is.null(pop)) res <- cbind.data.frame(pop=pop,res)
      return(res) # exit here
  }

  ## codom case ##
  # make one table by locus from x@tab
  kX <- seploc(x,res.type="matrix")

  if (oneColPerAll & all(x@ploidy == x@ploidy[1])){
    sep <- "/"
  }
  ## function to recode a genotype in form "A1[sep]...[sep]Ak" from frequencies
  recod <- function(vec,lab){
      if(any(is.na(vec))) return(NA)
      res <- paste(rep(lab,vec), collapse=sep)
      return(res)
  }


  # kGen is a list of nloc vectors of genotypes
  kGen <- lapply(1:length(kX), function(i) apply(kX[[i]],1,recod,x@all.names[[i]]))
  names(kGen) <- locNames(x)

  ## if use one column per allele
  if(oneColPerAll){
    if (all(x@ploidy == x@ploidy[1])){
      ## Build each row from allele counts: do not join and re-split labels.
      res <- lapply(seq_along(kX), function(i) {
          counts <- kX[[i]]
          out <- matrix(NA_character_, nrow=nrow(counts), ncol=x@ploidy[1])
          for (j in seq_len(nrow(counts))) {
              copies <- counts[j, ]
              if (anyNA(copies)) next
              if (sum(copies) > x@ploidy[1])
                  stop("Allele counts exceed ploidy for individual ",
                       indNames(x)[j], " at locus ", locNames(x)[i])
              observed <- rep(x@all.names[[i]], copies)
              out[j, seq_along(observed)] <- observed
          }
          out
      })
      res <- data.frame(res,stringsAsFactors=FALSE)
      names(res) <- paste(rep(locNames(x),each=x@ploidy[1]), 1:x@ploidy[1], sep=".")

      ## handle pop here
      if(!is.null(pop) & usepop) res <- cbind.data.frame(pop,res,stringsAsFactors=FALSE)
      rownames(res) <- indNames(x)
      return(res) # exit here
    } else {
      warning("All ploidies must be equal in order to separate the alleles.\nReturning one column per locus")
    }
  } # end if oneColPerAll

  ## build the final data.frame
  ## res <- as.data.frame(cbind(kGen,stringsAsFactors=FALSE))
  ## faster option
  res <- as.data.frame(do.call(cbind, kGen), stringsAsFactors=FALSE)
  rownames(res) <- indNames(x)
  colnames(res) <- locNames(x)

  ## handle pop here
  if(!is.null(pop) & usepop) res <- cbind.data.frame(pop,res,stringsAsFactors=FALSE)

  return(res)
}
