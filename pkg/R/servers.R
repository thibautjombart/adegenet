

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
