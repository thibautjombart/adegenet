

#' Web servers for adegenet
#'
#' The function \code{adegenetServer} opens up a web page providing a simple
#' user interface for some of the functionalities implemented in adegenet.
#' These servers have been developed using the package \code{shiny}.\cr
#'
#' Currently available servers include: \itemize{ \item \code{DAPC}: a server
#' for the Discriminant Analysis of Principal Components (see ?dapc) }
#'
#'
#' @aliases adegenetServer
#' @param what a character string indicating which server to start; currently
#' accepted values are: "DAPC"
#' @return The function invisibly returns NULL.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk} Caitlin Collins
#' @seealso \link{dapc}
#' @examples
#'
#' \dontrun{
#' ## this opens a web page for DAPC
#' adegenetServer()
#' }
#'

##################
## adegenetServer
##################
## function starting a web-server
## for adegenet tools
adegenetServer <- function(what=c("DAPC")){
    what <- match.arg(what)
    if(what=="DAPC"){
        .dapcServer()
    }

    return(invisible())
}


###############
## .dapcServer
###############
## hidden function - DAPC server
.dapcServer <- function(){
    runApp(system.file("dapcServer",package="adegenet"))
}




#######################
## .render.server.info
#######################
## INVISIBLE FUNCTION RENDERING SERVER INFO ##
.render.server.info <- function(){
    renderPrint(
            {
                cat("\n== R version ==\n")
                print(R.version)

                cat("\n== Date ==\n")
                print(date())

                cat("\n== adegenet version ==\n")
                print(packageDescription("adegenet", fields=c("Package", "Version", "Date", "Built")))

                cat("\n== shiny version ==\n")
                print(packageDescription("adegenet", fields=c("Package", "Version", "Date", "Built")))

                cat("\n== attached packages ==\n")
                print(search())
            }
                ) # end renderPrint
} # end .render.server.info
