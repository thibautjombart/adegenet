## First.lib <- function (lib, pkg){
## #.initAdegenetClasses()
## #.initAdegenetUtils()
##     library.dynam("adegenet", pkg, lib)
##     pkg.version <- packageDescription("adegenet", fields = "Version")

##     startup.txt <- paste("   ==========================\n    adegenet", pkg.version, "is loaded\n   ==========================\n\n - to start, type '?adegenet'\n - to browse adegenet website, type 'adegenetWeb()'\n - to post questions/comments: adegenet-forum@lists.r-forge.r-project.org\n\n")

##     packageStartupMessage(startup.txt)
## }

.onAttach <- function(libname, pkgname){
    pkg.version <- packageDescription("adegenet", fields = "Version")

    startup.txt <- paste("   ==========================\n    adegenet", pkg.version, "is loaded\n   ==========================\n\n - to start, type '?adegenet'\n - to browse adegenet website, type 'adegenetWeb()'\n - to post questions/comments: adegenet-forum@lists.r-forge.r-project.org\n\n")

    packageStartupMessage(startup.txt)
}
