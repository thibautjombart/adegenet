.onAttach <- function(libname, pkgname){
    # adegenet specific options -----------------------------------------------
    op <- options()
    op.adegenet <- list(
      adegenet.testcon = stdin() # for readLines, read from stdin. This allows it to be changed for tests.
    )
    toset <- !(names(op.adegenet) %in% names(op))
    if(any(toset)) options(op.adegenet[toset])

    # startup message ---------------------------------------------------------
    pkg.version <- packageDescription("adegenet", fields = "Version")
    startup.txt <- paste("\n   /// adegenet ", pkg.version, " is loaded ////////////",
                         "\n\n   > overview: '?adegenet'",
                         "\n   > tutorials/doc/questions: 'adegenetWeb()' ",
                         "\n   > bug reports/feature requests: adegenetIssues()\n\n", sep="")
    packageStartupMessage(startup.txt)
}
