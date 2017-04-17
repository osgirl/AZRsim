###############################################################################
###############################################################################
# This file contains general auxiliary functions
###############################################################################
###############################################################################

###############################################################################
# as_numeric: Converts a vector of elements to numeric elements
###############################################################################
#' Converts a vector of elements to numeric elements
#'
#' @param input A vector
#' @return A vector where the numeric elements are numeric and the non-numeric ones are NA
#' @examples
#' as_numeric(c(1,2,3))
#' as_numeric(c('1',2,3))
#' as_numeric(c('1','a',3))
#' as_numeric(c('1','2',NA))
#' as_numeric(factor(c(1, "4", 3)))
#' @export

as_numeric <- function(input) {
  suppressWarnings(as.numeric(as.character(input)))
}

###############################################################################
# isnumericVector: checks if a vector has only numeric entries
#                 "2" counts as numeric, "a" does not. "NA" does NOT count as
#                 numeric!
###############################################################################
#' Checks if a vector consists of numeric entries
#'
#' Checks if a vector has only numeric entries "2" counts as numeric, "a" does
#' not. "NA" and NA does NOT count as numeric either.
#'
#' @param input A vector
#' @return TRUE if numeric vector, FALSE if not
#' @examples
#' isnumericVector(c(1,2,3))
#' isnumericVector(c('1',2,3))
#' isnumericVector(c('1','a',3))
#' isnumericVector(c('1','2',NA))
#' @export

isnumericVector <- function(input) {
  return (!(NA %in% as_numeric(input)))
}


