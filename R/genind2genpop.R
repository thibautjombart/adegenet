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
#' - "0": if all genotypes of a population for a given allele are missing, count
#' value will be 0\cr
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
#' @param pop a factor giving the population of each genotype in 'x' OR a
#'   formula specifying which strata are to be used when converting to a genpop
#'   object. If none provided, population factors are sought in x@@pop, but if
#'   given, the argument prevails on x@@pop.
#' @param quiet logical stating whether a conversion message must be printed 
#'   (TRUE,default) or not (FALSE).
#' @param process.other a logical indicating whether the \code{@@other} slot 
#'   should be processed (see details).
#' @param other.action a function to be used when processing the \code{@@other} 
#'   slot. By default, 'mean' is used.
#' @return A genpop object. The component @@other in 'x' is passed to the 
#'   created genpop object.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \linkS4class{genind}, \linkS4class{genpop}
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
#' \dontrun{
#' data(microbov)
#' strata(microbov) <- data.frame(other(microbov))
#' summary(genind2genpop(microbov)) # Conversion based on population factor
#' summary(genind2genpop(microbov, ~coun)) # Conversion based on country
#' summary(genind2genpop(microbov, ~coun/spe)) # Conversion based on country and species
#' 
#' }
#'
#' @export genind2genpop
genind2genpop <- function(x, pop = NULL, quiet = FALSE, 
                          process.other = FALSE, other.action = mean){

    ## CHECKS ##
    if (!is.genind(x)) stop("x is not a valid genind object")
    checkType(x)
    if (!all(ploidy(x)[1]==ploidy(x))) stop("conversion to genpop not supported for varying ploidy")
    if (!is.null(pop)){
      if (is.language(pop)){
        setPop(x) <- pop
      } else {
        pop(x) <- pop
      }
    }

    if (is.null(pop(x))) {
        if(!quiet) warning("\npop is not provided either in x or in pop - assuming one single group")
        pop(x) <- factor(rep(1, nInd(x)))
    }

    if (!quiet) cat("\n Converting data from a genind to a genpop object... \n")

    ## tabcount is a matrix pop x alleles, counting alleles per pop
    tabcount <- apply(tab(x), 2, tapply, pop(x), sum, na.rm=TRUE)

    ## restitute matrix class when only one pop
    if(is.null(dim(tabcount))) {
        lab.col <- names(tabcount)
        tabcount <- matrix(tabcount,nrow=1)
        colnames(tabcount) <- lab.col
        rownames(tabcount) <- levels(pop(x))[1]
    }

    ## MAKE FINAL OBJECT ##
    prevcall <- match.call()
    res <- new("genpop", tab=tabcount, prevcall=prevcall, ploidy=x@ploidy[1], type=x@type)

    ## handle @other here
    res@other <- x@other
    if(process.other){
        ## auxiliary function doing the job
        fOther <- function(e){
            N <- nrow(x@tab)
            if(is.vector(e) && is.numeric(e) && length(e)==N){ # numeric vector
                res <- tapply(e, pop(x), other.action)
                return(res)
            } else if(is.matrix(e) && is.numeric(e) && nrow(e)==N){ # numeric matrix
                res <- apply(e, 2, function(vec) tapply(vec, pop(x), other.action))
                colnames(res) <- colnames(e)
                return(res)
            } else if(is.data.frame(e) && nrow(e)==N && all(sapply(e,is.numeric)) ){ # df of numeric vectors
                res <- lapply(e, function(vec) tapply(vec, pop(x), other.action))
                res <- data.frame(res)
                names(res) <- names(e)
                return(res)
            } else return(e)
        } # end fOther

        res@other <- lapply(res@other, fOther)
    } # end if(process.other)

    if(!quiet) cat("\n...done.\n\n")

    return(res)

} # end genind2genpop
