#' @title Distance matrices from genlight objects
#' 
#' @name dist.genlight
#' @rdname dist.genlight
#' 
#' @description 
#' Create distance matrices from genlight objects
#' 
#' @details 
#' 
#' The creation of distance matrices, matrices of numbers that describe how different each of the samples are, is a fundamental task in the statistical analysis of individuals or populations (i.e., groups of individuals).
#' However, there isn't actually a function that creates distance matrices from genlight objects in adegenet.
#' Instead, the authors of adegenet created an `as.matrix()` function that converts a genlight object to a matrix.
#' This is clever because the function `dist()` in the package `stats` tries to convert whatever object it is given to a matrix.
#' The result is that when you call `dist()` on a genlight object it uses the `dist()` function to create a distance matrix.
#' The reason this is clever is because it uses pre-existing code.
#' The downside is that because there is no function to specifically create distance matrices from genlight objects in adegenet, there is no documentation in genlight for how this is done.
#' And because the author of `dist()` never anticipated it could be used on genlight objects, there is no documentation for it there either.
#' And we can find documentation for this function with `?dist`.
#' To summarize, we can create a distance matrix from a genlight object using `dist()`.
#' 
#'  
#' There are also functions to create distance matrices from genlight objects that exist in other packages.
#' The function `bitwise.dist()` in the package [poppr](https://CRAN.R-project.org/package=poppr) is an example.
#' We can find documentation for this function with `?poppr::bitwise.dist`.
#' Again, a downside of this is that you need to know where to look for this information or you may not find it.
#' 
#' 
#' Lastly, because you can use `as.matrix()` on your genlight object, and most distance algorithms can use this matrix as input, you can use this as an intermediate step to create a matrix from your genlight object and pass it to your distance algorithm of choice.
#' Options include [ade4](https://CRAN.R-project.org/package=ade4), `vegdist()` in [vegan](https://CRAN.R-project.org/package=vegan), or `daisy()` in [cluster](https://CRAN.R-project.org/package=cluster).
#' Note that it is up to you to determine which distance metric is best for your analysis.
#' A number of options therefore exist for creating distance matrices from genlight objects.

