###############################################################################
###############################################################################
# This file contains some AZRmodel math functions
###############################################################################
###############################################################################


###############################################################################
# gt
###############################################################################
# Greater than function
#
# Needed for AZRmodels especially in events and piecewise statements
#
# @param x First value
# @param y Second value
# @return as.numeric(x>y)
# @examples
# gt(4,3)
# @export

gt <- function (x,y) {
  return(as.numeric(x>y))
}

###############################################################################
# ge
###############################################################################
# Greater equal function
#
# Needed for AZRmodels especially in events and piecewise statements
#
# @param x First value
# @param y Second value
# @return as.numeric(x>=y)
# @examples
# ge(4,3)
# @export

ge <- function (x,y) {
  return(as.numeric(x>=y))
}

###############################################################################
# lt
###############################################################################
# Less than function
#
# Needed for AZRmodels especially in events and piecewise statements
#
# @param x First value
# @param y Second value
# @return as.numeric(x<y)
# @examples
# lt(4,3)
# @export

lt <- function (x,y) {
  return(as.numeric(x<y))
}

###############################################################################
# le
###############################################################################
# Less equal function
#
# Needed for AZRmodels especially in events and piecewise statements
#
# @param x First value
# @param y Second value
# @return as.numeric(x<=y)
# @examples
# le(4,3)
# @export

le <- function (x,y) {
  return(as.numeric(x<=y))
}

###############################################################################
# mod
###############################################################################
# Modulo/mod function
#
# @param x First value
# @param y Second value
# @return x %% y
# @examples
# mod(4,3)
# mod(7,2)
# @export

mod <- function (x,y) {
  return(x %% y)
}

###############################################################################
# and
###############################################################################
# and function
#
# Needed for AZRmodels especially in events and piecewise statements
#
# @param ... several numeric input arguments
# @return numeric "and" connection of input arguments
# @examples
# and(1)
# and(1,0)
# and(0,1)
# and(0,0)
# and(0,0,0)
# and(0,0,1)
# and(1,2,1)
# @export

and <- function (...) {
  varargin <- unlist(list(...))
  nargin <- length(varargin)
  result <- TRUE
  for (k in 1:nargin) {
    result <- result & varargin[k]
  }
  return(as.numeric(result))
}

###############################################################################
# or
###############################################################################
# or function
#
# Needed for AZRmodels especially in events and piecewise statements
#
# @param ... several numeric input arguments
# @return numeric "or" connection of input arguments
# @examples
# or(1)
# or(1,0)
# or(0,1)
# or(0,0)
# or(0,0,0)
# @export

or <- function (...) {
  varargin <- unlist(list(...))
  nargin <- length(varargin)
  result <- FALSE
  for (k in 1:nargin) {
    result <- result | varargin[k]
  }
  return(as.numeric(result))
}

###############################################################################
# multiply
###############################################################################
# multiply function (MathML)
#
# Needed for AZRmodels especially in some SBML import cases
#
# @param ... several numeric input arguments
# @return multiplication of all input arguments
# @examples
# multiply(3)
# multiply(1,2)
# multiply(0,1)
# multiply(5,4)
# multiply(1,2,3)
# multiply(3,2,1)
# multiply(1,2,1)
# @export

multiply <- function (...) {
  varargin <- unlist(list(...))
  nargin <- length(varargin)
  result <- 1
  for (k in 1:nargin) {
    result <- result * varargin[k]
  }
  return(result)
}

###############################################################################
# piecewise
###############################################################################
# Piecewise function
#
# This function implements support for the SBML / MATHML piecewise operator.
# It is additionally a convenient way for the implementation of IF ELSEIF ELSEIF
# ELSE statements in AZRmodels.
#
# result <- piecewise(resultiftrue1,decision1,resultiftrue2,decision2,...,resultiftruen,decisionn)
# result <- piecewise(resultiftrue1,decision1,resultiftrue2,decision2,...,resultiftruen,decisionn,defaultresult)
#
# @param ... Decision and result values
# @return Value
# @examples
# piecewise(3,lt(4,3),1)
# piecewise(3,lt(3,4),1)
# piecewise(3,ge(3,4),1,ge(2,3),0)
# piecewise(3,ge(5,4),1,ge(4,3),0)
# piecewise(3,ge(3,4),1,ge(4,3),0)
# piecewise(3,ge(3,4),1,ge(4,3))
# piecewise(3,ge(3,4),1,ge(2,3),3)
# @export

piecewise <- function (...) {
  varargin <- unlist(list(...))
  nargin <- length(varargin)
  result <- NULL

  # check if odd or even number of input arguments
  oddnumber = mod(nargin,2)

  for (k in seq(1,nargin-oddnumber,2)) {
    if (varargin[k+1] != 0) {
      result <- varargin[k]
      return(result)
    }
  }

  if (is.null(result)) {
    if (oddnumber==1) {
      result <- varargin[nargin]
    } else {
      stop('piecewise: piecewise statement is wrongly defined - missing (but needed) default value');
    }
  }

  return(result)
}

###############################################################################
# Interpolation functions
# interp0
# interp1
# interpcs
# interpcse
# Important: due to R syntax the use of these functions in AZRmodels is limited
# to C-code simulation. Therefor, when simulating a model with any of these
# functions using deSolve, only an error will be shown and thus the functions
# defined below are merely "dummy" functions returning an error.
###############################################################################

# interp0: zero order interpolation function
#
# zero order interpolation function (lookup table)
# If of limits then the extreme points in y are taken as output.
#
# @param x vector of function arguments
# @param y vector of function values at the points given by x
# @param xi scalar value for which to determine y by zero order interpolation
# @return Interpolated value
# @export
interp0 <- function (x,y,xi) {
  if (xi >= x[length(x)]) return(y[length(y)])
  if (xi <= x[1]) return(y[1])
  return(stats::approx(x,y,xi,method="const")$y)
}

# interp1: first order (linear) interpolation function
#
# linear interpolation function (lookup table)
# If of limits then the extreme points in y are taken as output.
#
# @param x vector of function arguments
# @param y vector of function values at the points given by x
# @param xi scalar value for which to determine y by linear interpolation
# @return Interpolated value
# @export
interp1 <- function (x,y,xi) {
  if (xi >= x[length(x)]) return(y[length(y)])
  if (xi <= x[1]) return(y[1])
  return(stats::approx(x,y,xi,method="linear")$y)
}

# interpcs: cubic spline interpolation function
#
# cubic spline interpolation function
# If of limits then the extreme points in y are taken as output.
# Note this function uses a slightly different approach as the corresponding
# C-code implementation, leading to slightly different results.
#
# @param x vector of function arguments
# @param y vector of function values at the points given by x
# @param xi scalar value for which to determine y by cubic spline interpolation
# @return Interpolated value
# @export
interpcs <- function (x,y,xi) {
  if (xi >= x[length(x)]) return(y[length(y)])
  if (xi <= x[1]) return(y[1])
  return(stats::spline(x=x,y=y,xout=xi,method="natural")$y)
}

# interpcse: Cubic spline interpolation with endpoints
# Implemented as C code function in the src folder.


