# 2021-06-05
#
# This script will re-build the adegenet tutorials on your machine and throw a
# message got STDERR if anything goes wrong with one of the tutorials. 
#
# If you don't have tinytex, it will be installed for you.
if (!require("tinytex")) {
  message("installing tinytex because you will need it")
  install.packages("tinytex", repos = "https://cloud.r-project.org")
  tinytex::install_tinytex()
}

files <- list.files(pattern = "Rnw")
if (length(files) == 0) {
  wd <- setwd(file.path(".", "tutorials"))
  on.exit(setwd(wd), add = TRUE)
  files <- list.files(pattern = "Rnw")
}

library("knitr")
for (f in files) {
  tryCatch(knit2pdf(f), 
    error = function(e) message("ERROR\n", e$message)
  )
}

