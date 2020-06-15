#' @import methods
#' @import stats
#' @import data.table
#' @useDynLib variantprobs, .registration = TRUE

NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("mypackage", libpath)
}
