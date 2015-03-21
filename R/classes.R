########################################################################
# adegenet classes definitions. All classes are S4.
#
# Thibaut Jombart, November 2007
# t.jombart@imperial.ac.uk
########################################################################

###############################
# Two classes of R object are
# defined :
# gen - common part to genind and genpop
# genind - genotypes of individuals
# genpop - allelic frequencies of populations
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


####################
# virtual class gen
####################
.gen.valid <- function(object){
  # this function tests only the consistency
  # of the length of each component
  p <- ncol(object@tab)
  k <- length(unique(object@loc.names))


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

  if(!is.null(object@loc.nall)){
      if(length(object@loc.nall) != k) {
          cat("\ninvalid length in loc.nall\n")
          return(FALSE)
      }
  }

  temp <- table(object@loc.names[object@loc.names!=""])
  if(any(temp>1)) {
      warning("\nduplicate names in loc.names:\n")
      print(temp[temp>1])
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
                               loc.names = "character",
                               loc.fac = "factorOrNULL",
                               loc.nall = "intOrNum",
                               all.names = "listOrNULL",
                               call = "callOrNULL",
                               "VIRTUAL"),
         prototype(tab=matrix(0L, ncol=0,nrow=0), loc.nall=integer(0), call=NULL))

setValidity("gen", .gen.valid)





########################
# virtual class indInfo
########################
setClass("indInfo", representation(ind.names = "character",
                                   pop = "factorOrNULL",
                                   pop.names = "charOrNULL",
                                   ploidy = "integer",
                                   type = "character",
                                   other = "listOrNULL", "VIRTUAL"),
         prototype(pop=NULL, pop.names = NULL, type = "codom", ploidy = as.integer(2), other = NULL))





###############
# Class genind
###############
.genind.valid <- function(object){
    if(!.gen.valid(object)) return(FALSE)

    if(length(object@ind.names) != nrow(object@tab)) {
        cat("\ninvalid length in ind.names\n")
        return(FALSE)
    }

    temp <- table(object@ind.names[object@ind.names!=""])
    if(any(temp>1)) {
        warning("\nduplicate names in ind.names:\n")
        print(temp[temp>1])
    }

    if(typeof(object@tab)!="integer"){
        warning("@tab does not contain integers; as of adegenet_1.5-0, numeric values are no longer used")
        ## cat("\ntab does not contain integers; as of adegenet_1.5-0, numeric values are no longer used")
        ## return(FALSE)
    }


    if(!is.null(object@pop)){ # check pop

        if(length(object@pop) != nrow(object@tab)) {
            cat("\npop is given but has invalid length\n")
            return(FALSE)
        }

        if(is.null(object@pop.names)) {
            cat("\npop is provided without pop.names")
        }


        if(length(object@pop.names) != length(levels(object@pop))) {
            cat("\npop.names has invalid length\n")
            return(FALSE)
        }

        temp <- table(object@pop.names[object@pop.names!=""])
        if(any(temp>1)) {
            warning("\nduplicate names in pop.names:\n")
            print(temp[temp>1])
        }

    } # end check pop

    # Check population hierarchy
    if (!is.null(object@hierarchy)){
      if (nrow(object@hierarchy) != nrow(object@tab)){
        cat("\na hierarchy is defined has invalid length\n")
        return(FALSE)
      }

      dups <- duplicated(colnames(object@hierarchy))
      if (any(dups)){
        cat("\nduplicated names found in @hierarchy slot:\n")
        dups <- colnames(object@hierarchy)[dups]
        cat(paste0(dups, collapse = ", "))
        return(FALSE)
      }
    }

    ## check ploidy
    if(any(object@ploidy < 1L)){
        cat("\nploidy inferior to 1\n")
        return(FALSE)
    }
    if(length(object@ploidy)!=nInd(object)){
        warning("as of adegenet_1.5-0, @ploidy should contain one value per individual")
    }

    ## check type of marker
    if(!object@type %in% c("codom","PA") ){
        cat("\nunknown type of marker\n")
        return(FALSE)
    }


    return(TRUE)
} #end .genind.valid

setClass("genind", contains=c("gen", "indInfo"), 
          representation = representation(hierarchy = "dfOrNULL"))
setValidity("genind", .genind.valid)



########################
# virtual class popInfo
########################
setClass("popInfo", representation(pop.names = "character", ploidy = "integer",
                                   type = "character", other = "listOrNULL", "VIRTUAL"),
         prototype(type = "codom", ploidy = as.integer(2), other = NULL))



###############
# Class genpop
###############
.genpop.valid <- function(object){
    if(!.gen.valid(object)) return(FALSE)
    if(length(object@pop.names) != nrow(object@tab)) {
        cat("\ninvalid length in pop.names\n")
        return(FALSE)
    }

    temp <- table(object@pop.names[object@pop.names!=""])
    if(any(temp>1)) {
        warning("\nduplicate names in pop.names:\n")
        print(temp[temp>1])
    }

     ## check ploidy
    if(object@ploidy < 1L){
        cat("\nploidy inferior to 1\n")
        return(FALSE)
    }

    ## check type of marker
    if(!object@type %in% c("codom","PA") ){
        cat("\nunknown type of marker\n")
        return(FALSE)
    }

    return(TRUE)
} #end .genpop.valid

setClass("genpop", contains=c("gen", "popInfo"))
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

