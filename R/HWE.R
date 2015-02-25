##################
# HWE.test.genind
##################



#' Hardy-Weinberg Equilibrium test for multilocus data
#' 
#' The function \code{HWE.test} is a generic function to perform Hardy-Weinberg
#' Equilibrium tests defined by the \code{genetics} package. adegenet proposes
#' a method for \code{genind} objects.\cr
#' 
#' The output can be of two forms:\cr - a list of tests (class \code{htest})
#' for each locus-population combinaison \cr - a population x locus matrix
#' containing p-values of the tests
#' 
#' Monte Carlo procedure is quiet computer-intensive when large datasets are
#' involved. For more precision on the performed test, read \code{HWE.test}
#' documentation (\code{genetics} package).
#' 
#' @param x an object of class \code{genind}.
#' @param pop a factor giving the population of each individual. If NULL, pop
#' is seeked from x\$pop.
#' @param permut a logical passed to \code{HWE.test} stating whether Monte
#' Carlo version (TRUE) should be used or not (FALSE, default).
#' @param nsim number of simulations if Monte Carlo is used (passed to
#' \code{HWE.test}).
#' @param hide.NA a logical stating whether non-tested loci (e.g., when an
#' allele is fixed) should be hidden in the results (TRUE, default) or not
#' (FALSE).
#' @param res.type a character or a character vector whose only first argument
#' is considered giving the type of result to display. If "full", then a list
#' of complete tests is returned. If "matrix", then a matrix of p-values is
#' returned.
#' @return Returns either a list of tests or a matrix of p-values. In the first
#' case, each test is designated by locus first and then by population. For
#' instance if \code{res} is the "full" output of the function, then the test
#' for population "PopA" at locus "Myloc" is given by res$Myloc$PopA. If
#' \code{res} is a matrix of p-values, populations are in rows and loci in
#' columns. P-values are given for the upper-tail: they correspond to the
#' probability that an oberved chi-square statistic as high as or higher than
#' the one observed occured under H0 (HWE).\cr
#' 
#' In all cases, NA values are likely to appear in fixed loci, or entirely
#' non-typed loci.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link[genetics]{HWE.test}},\code{\link[stats]{chisq.test}}
#' @keywords manip multivariate
#' @examples
#' 
#' \dontrun{
#' data(nancycats)
#' obj <- nancycats
#' if(require(genetics)){
#' obj.test <- HWE.test(obj)
#' 
#' # pvalues matrix to have a preview
#' HWE.test(obj,res.type="matrix")
#' 
#' #more precise view to...
#' obj.test$fca90$P10
#' }
#' }
#' 
#' @export HWE.test.genind
HWE.test.genind <- function(x,pop=NULL,permut=FALSE,nsim=1999,hide.NA=TRUE,res.type=c("full","matrix")){

  if(!is.genind(x)) stop("x is not a valid genind object")
  if(x@ploidy != as.integer(2)) stop("not implemented for non-diploid genotypes")
  checkType(x)

  if(!require(genetics)) stop("genetics package is required. Please install it.")
  if(is.null(pop)) pop <- x@pop
  if(is.null(pop)) pop <- as.factor(rep("P1",nrow(x@tab)))
  res.type <- tolower(res.type[1])
  if(res.type != "full" && res.type != "matrix") stop("unknown res.type specified.")

  kGen <- genind2genotype(x,pop=pop,res.type="list")

  # ftest tests HWE for a locus and a population
  ftest <- function(vec,permut=permut,nperm=nsim){
    temp <- unique(vec)
    temp <- temp[!is.na(temp)]
    if(length(temp) < 2) return(NA)
    if(res.type=="full") {
      res <- HWE.chisq(vec, simulate.p.value=permut, B=nperm)
    } else {
      res <- HWE.chisq(genotype(vec), simulate.p.value=permut, B=nperm)$p.value
    }
    return(res)
  }

  res <- lapply(kGen,function(e) lapply(e,ftest,permut,nsim))

  # clean non-tested elements in the results list
  if(hide.NA && res.type=="full"){
    newres=list()
    tokeep <- which(unlist(lapply(res,function(e) !all(is.na(e)))))
    if(length(tokeep) > 0) for(i in 1:length(tokeep)) {newres[[i]] <- res[[tokeep[i]]]}
    newres <- lapply(newres,function(e) {e[!is.na(e)] })
    names(newres) <- names(res)[tokeep]
    res <- newres
  }

  if(res.type=="matrix"){
    res <- as.data.frame(lapply(res,unlist))
    rnam <- rownames(res)
    rownames(res) <- gsub(".X-squared","",rnam)
    res <- as.matrix(res)
  }

  return(res)
}
