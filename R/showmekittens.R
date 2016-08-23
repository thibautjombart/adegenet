#' When you need a break...
#'
#' Genetic data analysis can be a harsh, tiring, daunting task.
#' Sometimes, a mere break will not cut it.
#' Sometimes, you need a kitten.
#'
#' \author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#'
#' @details
#'
#' Please send us more! Either pull request or submit an issue with a URL (use
#' \code{adegenetIssues()}).
#'
#'
#' @param x the name or index of the video to display; if NULL, a random video is chosen
#'
#' @param list a logical indicating if the list of available videos should be displayed
#'
showmekittens <- function(x = NULL, list = FALSE){
    ## 'pool' is a named character vector of video URLs
    pool <- c(capucine = "http://www.youtube.com/watch?v=KIePsbJSS04",
              vacuum = "https://www.youtube.com/watch?v=uiyKVWqxXWM")

    ## either we return the list of videos, or we show one
    if (list) {
        return(pool)
    }

    if (is.null(x)) {
        x <- sample(seq_along(pool), 1L)
    }

    ## check that x is okay
    if (is.numeric(x) && (x < 1 || x > length(pool))) {
        stop(sprintf("Video index (%d) is wrong; there are currently %d videos in the list", x, length(pool)))
    }

    if (is.character(x) && !x %in% names(pool)) {
        stop(sprintf("Video name (%s) is not in the list; use the option 'list=TRUE' to see available videos.", x))
    }

    browseURL(pool[x])
}
