#' Print function for FEmrt
#'
#' Print the results of a FEmrt object
#'
#' @param x fitted tree of class \code{FEmrt}.
#' @param \dots additional arguments to be passed.
#' @details
#' The function returns the objects concerning the analysis results.
#'
#' @examples data(SimData)
#' test <- FEmrt(efk~m1+m2+m3+m4+m5, vark, data=SimData, c=1)
#' print(test)
#' @export
print.FEmrt<- function(x, ...){
  if (length(x$n) < 2) {
    cat("\n")
    cat("Fixed Effects Meta-Tree (K = ", sum(x$n), " studies); ",
        sep = "")
    cat("\n")
    print(x$call)
    cat("\n")
    cat("No moderator effect was detected" )


  } else {
    cat("\n")
    cat("Fixed Effects Meta-tree (K = ", sum(x$n), " studies); ",
        sep = "")
    cat("\n")
    print(x$call)
    cat("\n")
    cat("A tree with ", length(x$n), " terminal nodes was detected", sep="" )
    cat("\n")
    cat("The moderators are ", paste(as.character(x$moderators), collapse = ", "), sep = "")
  }
}
