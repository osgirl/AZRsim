###############################################################################
###############################################################################
# This file contains general auxiliary functions
###############################################################################
###############################################################################

###############################################################################
# veclocate: returns indices of elements in a vector that match a provided criterion
###############################################################################
# Returns indices of elements in a vector that match a provided criterion
#
# @param input A criterion applied to a vector
# @return A vector with indices of matching elements
# @examples
# veclocate(c(1,4,5)>=2)  # should return c(2,3)
# veclocate(c(1,4,5)==3)  # should return numeric(0)
# veclocate(c("abc",4,5)=="abc")  # should return c(1)
# @export

# veclocate <- function(input) {
#   x <- as.numeric(input)
#   y <- x*c(1:length(x))
#   return(y[y!=0])
# }

veclocate <- function(input) {
  which(input)
}

###############################################################################
# isnumericVector: checks if a vector has only numeric entries
#                 "2" counts as numeric, "a" does not. "NA" does NOT count as
#                 numeric!
###############################################################################
# Checks if a vector consists of numeric entries
#
# Checks if a vector has only numeric entries "2" counts as numeric, "a" does
# not. "NA" and NA does NOT count as numeric either.
#
# @param input A vector
# @return TRUE if numeric vector, FALSE if not
# @examples
# isnumericVector(c(1,2,3))  # should return TRUE
# isnumericVector(c('1',2,3))  # should return TRUE
# isnumericVector(c('1','a',3))  # should return FALSE
# isnumericVector(c('1','2',NA))  # should return FALSE

# isnumericVector <- function(input) {
#   if (sum(as.numeric(is.na(deString(input)))) > 0)
#     return(FALSE)
#   else
#     return(TRUE)
# }

isnumericVector <- function(input) {
  input <- suppressWarnings(as.numeric(input))
  return(!any(is.na(input)))
}
