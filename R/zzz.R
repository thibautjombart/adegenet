.onAttach <- function(libname, pkgname){
    pkg.version <- packageDescription("adegenet", fields = "Version")

    startup.txt <- paste("\n   /// adegenet ", pkg.version, " is loaded ////////////",
                         "\n\n   > overview: '?adegenet'",
                         "\n   > tutorials/doc/questions: 'adegenetWeb()' ",
                         "\n   > bug reports/feature requests: adegenetIssues()\n\n", sep="")

    packageStartupMessage(startup.txt)
}
