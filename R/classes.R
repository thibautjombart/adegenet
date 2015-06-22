########################################################################
## adegenet classes definitions. All classes are S4.
##
## Initial development: Thibaut Jombart, November 2007
##
## Major reform for adegenet 2.0.0 (March-August 2015)
##
## t.jombart@imperial.ac.uk
########################################################################

###############################
# Two classes of R object are
# defined :
# gen - common part to genind and genpop
# genind - allele counts for individuals
# genpop - allele counts for populations
###############################


###############################################################
###############################################################
# AUXILIARY FUNCTIONS
###############################################################
###############################################################




###############################################################
###############################################################
# CLASSES DEFINITION
###############################################################
###############################################################

##.initAdegenetClasses <- function(){


####################
# Unions of classes
####################
setClassUnion("listOrNULL", c("list","NULL"))
setClassUnion("factorOrNULL", c("factor","NULL"))
setClassUnion("charOrNULL", c("character","NULL"))
setClassUnion("callOrNULL", c("call","NULL"))
setClassUnion("intOrNum", c("integer","numeric","NULL"))
setClassUnion("intOrNULL", c("integer","NULL"))
setClassUnion("dfOrNULL", c("data.frame", "NULL"))
setClassUnion("formOrNULL", c("formula", "NULL"))


####################
# virtual class gen
####################
.gen.valid <- function(object){
  # this function tests only the consistency
  # of the length of each component
  p <- ncol(object@tab)
  k <- length(levels(object@loc.fac))


  if(!is.null(object@loc.fac)){
      if(length(object@loc.fac) != p) {
          cat("\ninvalid length for loc.fac\n")
          return(FALSE)
      }

      if(length(levels(object@loc.fac)) != k) {
          cat("\ninvalid number of levels in loc.fac\n")
          return(FALSE)
      }
  }

  if(!is.null(object@loc.n.all)){
      if(length(object@loc.n.all) != k) {
          cat("\ninvalid length in loc.n.all\n")
          return(FALSE)
      }
  }

  if(!is.null(object@all.names)){
      if(length(unlist(object@all.names)) != p) {
          cat("\ninvalid length in all.names\n")
          return(FALSE)
      }
  }

  return(TRUE)

}# end .gen.valid


setClass("gen", representation(tab = "matrix",
                               loc.fac = "factorOrNULL",
                               loc.n.all = "intOrNum",
                               all.names = "listOrNULL",
                               ploidy = "integer",
                               type = "character",
                               other = "listOrNULL",
                               call = "callOrNULL",
                               "VIRTUAL"),
         prototype(tab=matrix(0L, ncol=0,nrow=0),
                   loc.fac=NULL,
                   loc.n.all=integer(0),
                   all.names=NULL,
                   ploidy=integer(0),
                   type=character(0),
                   other=NULL,
                   call=NULL))

setValidity("gen", .gen.valid)




########################
# virtual class indInfo
########################
setClass("indInfo", representation(pop = "factorOrNULL",
                                   strata = "dfOrNULL",
                                   hierarchy = "formOrNULL",
                                   "VIRTUAL"),
         prototype(pop=NULL, stata=NULL, hierarchy=NULL))





###############
# Class genind
###############

setClass("genind", contains=c("gen", "indInfo"))

.genind.valid <- function(object){

    validation <- TRUE
    if(!.gen.valid(object)) return(FALSE)

    if(typeof(object@tab)!="integer"){
        warning("@tab does not contain integers; as of adegenet_2.0-0, numeric values are no longer used")
        ## message("\ntab does not contain integers; as of adegenet_1.5-0, numeric values are no longer used")
        ## validation <- FALSE
    }


    if(!is.null(object@pop)){ # check pop

        if(length(object@pop) != nrow(object@tab)) {
            message("\npop is given but has invalid length\n")
            validation <- FALSE
        }

    } # end check pop

    # Check population strata
    if (!is.null(object@strata)){
      if (nrow(object@strata) != nrow(object@tab)){
        message("\na strata is defined has invalid length\n")
        validation <- FALSE
      }

      dups <- duplicated(colnames(object@strata))
      if (any(dups)){
        message("\nduplicated names found in @strata slot:\n")
        dups <- colnames(object@strata)[dups]
        message(paste0(dups, collapse = ", "))
        validation <- FALSE
      }
    }

    # TODO: CHECK HIERARCHY FORMULA

    ## check ploidy
    if(any(object@ploidy < 1L)){
        message("\nploidy inferior to 1\n")
        validation <- FALSE
    }
    if(length(object@ploidy)!= nrow(object@tab)){
        warning("as of adegenet_2.0-0, @ploidy should contain one value per individual")
    }

    ## check type of marker
    if(!object@type %in% c("codom","PA") ){
        message("\nunknown type of marker\n")
        validation <- FALSE
    }


    return(validation)
} #end .genind.valid

setValidity("genind", .genind.valid)



########################
# virtual class popInfo
########################
setClass("genpop", contains=c("gen"))



###############
# Class genpop
###############
.genpop.valid <- function(object){

    validation <- TRUE

    if(!.gen.valid(object)) return(FALSE)

     ## check ploidy
    if(length(object@ploidy) > 1 && object@ploidy < 1L){
        message("\nploidy inferior to 1\n")
        validation <- FALSE
    }

    ## check type of marker
    if(!object@type %in% c("codom","PA") ){
        message("\nunknown type of marker\n")
        validation <- FALSE
    }

    return(validation)
} #end .genpop.valid

setValidity("genpop", .genpop.valid)







###############################################################
###############################################################
## MISCELLANEOUS METHODS
###############################################################
###############################################################



#################
# Function names
#################
setMethod("names", signature(x = "genind"), function(x){
    return(slotNames(x))
})

setMethod("names", signature(x = "genpop"), function(x){
    return(slotNames(x))
})

