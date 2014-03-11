

##################
## adegenetServer
##################
## function starting a web-server
## for adegenet tools
adegenetServer <- function(what=c("DAPC")){
    what <- tolower(what)
    what <- match.arg(what)
    if(what=="dapc"){
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
