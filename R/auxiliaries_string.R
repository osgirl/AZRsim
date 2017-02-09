###############################################################################
###############################################################################
# This file contains functions related to string manipulations
###############################################################################
###############################################################################


###############################################################
# strtrimM: truncates strings by removing leading and trailing
#             white spaces (shadows the base strtrimM function but
#             is more in line with strtrimM in PHP and MATLAB)
###############################################################
#' Truncates strings by removing leading and trailing white spaces
#'
#' Whitespace characters are the following: intToUtf8(c(9, 10, 11, 12, 13, 32))
#'
#' @param input A character string
#' @return Character string with leading and trailing whitespaces removed
#' @examples
#' strtrimM("abcde")
#' strtrimM("a b c d e")
#' strtrimM(" a b c d e    ")
#' strtrimM(intToUtf8(c(9, 10, 11, 12, 13, 32, 97,  32,  98,  32,  99,  32, 100,  32, 101, 32, 13)))
#' @export

strtrimM <- function(input) {
  return(gsub("^\\s+|\\s+$", "", input))
}


###############################################################
# strrepM: string replace function
###############################################################
#' String replacement function
#'
#' @param origstr A character string in which to replace sub strings
#' @param oldsubstr Search string to replace
#' @param newsubstr A string to replace the search string if found
#' @return Resulting output string
#' @examples
#' strrepM("Hello World!","World","Kitty")
#' strrepM("a b c d e"," ","")
#' @export

strrepM <- function(origstr,oldsubstr,newsubstr) {
  return(gsub(oldsubstr, newsubstr, origstr, fixed="TRUE"))
}


###############################################################
# strremWhite: Removes all white spaces in strings
###############################################################
#' Removes all white spaces in a string
#'
#' Whitespace characters are the following: intToUtf8(c(9, 10, 11, 12, 13, 32))
#'
#' @param inputText A character string
#' @return Character string with whitespaces removed
#' @examples
#' strremWhite("abcde")
#' strremWhite("a b c d e")
#' strremWhite(" a b c d e    ")
#' strremWhite(intToUtf8(c(9, 10, 11, 12, 13, 32, 97,  32,  98,  32,  99,  32, 100,  32, 101, 32, 13)))
#' @export

strremWhite <- function(inputText) {
  return(gsub("\\s+|\\s+", "", inputText))
}


###############################################################
# strmatch: Returns indices of elements in a vector of that match a search string
###############################################################
#' Returns indices of elements in a vector of that match a search string
#'
#' Exact, case sensitive matching is used.
#' Function can also be used for other objects than strings.
#'
#' @param inputVec A vector of strings
#' @param searchString A string to search for (case sensitive,
#'   exact matching) in inputVec
#' @return A vector with indices of elements that have an exact match.
#'   NULL is returned if string not found
#' @examples
#' strmatch(c("hello","test","a","e","hello"),"hello")
#' strmatch(c("hello","test","a","e","hello"),"Hello")
#' strmatch(c("hello","test","a","e","hello"),"xyz")
#' strmatch(c(1,2,3,4,8,3),3)
#' @export

strmatch <- function(inputVec,searchString) {
  x <- which(inputVec==searchString)
  if (length(x)==0) return (NULL)
  return(unname(x))
}


###############################################################
# strexplode: Splits strings by separator
###############################################################
#' Splits a string based on a separator substring
#'
#' By default a "," is used as separator.
#'
#' @param input A string to split
#' @param separator A substring to use as separator
#' @return A character vector with the split elements
#' @examples
#' strexplode(c("Hello,,Test,,X,,Z"))
#' strexplode(c("Hello,,Test,,X,,Z"),",")
#' strexplode(c("Hello,,Test,,X,,Z"),",,")
#' strexplode(c("Hello,,Test,,X,,Z"),"§")
#' @export

strexplode <- function(input,separator=",") {
  return(unlist(strsplit(input,separator)))
}

###############################################################
# strexplodePC: Splits strings by separator, ignore separator in parentheses
###############################################################
#' Splits strings by separator only if separators not in parentheses
#'
#' By default a "," is used as separator.
#'
#' @param input A string to split
#' @param separator A substring to use as separator
#' @param group "round", "square", or "curly" allowing grouping by parantheses
#'   of defined type in which elements are not exploded
#' @return A character vector with the split elements
#' @examples
#' strexplodePC("(a,b),(c,d)")
#' strexplodePC("((a,b),(c,d))") # compare with previous - lowest level considered
#' strexplodePC("(a,b),(c,d))")
#' strexplodePC("(a;b);(c;d)",";")
#' @export

strexplodePC <- function(input,separator=",",group="round") {
  if (group=="round") {
    groupStart <- "("
    groupEnd   <- ")"
  }
  if (group=="square") {
    groupStart <- "["
    groupEnd   <- "]"
  }
  if (group=="curly") {
    groupStart <- "{"
    groupEnd   <- "}"
  }
  if (group!="round" && group!="square" && group!="curly") stop("strexplodePC: wrong group definition")

  # DO THE EXPLOSION
  elements        <- c()
  openParenthesis <- 0
  lastIndex       <- 1
  elementIndex    <- 1

  for (k2 in 1:nchar(input)) {
    if (substr(input,k2,k2) == groupStart) {
      openParenthesis <- openParenthesis + 1
    } else {
      if (substr(input,k2,k2) == groupEnd) {
        openParenthesis <- openParenthesis - 1;
      } else {
        if ((substr(input,k2,k2) == separator) && (openParenthesis == 0)) {
          elements[elementIndex] <- strtrimM(substr(input,lastIndex,k2-1))
          elementIndex           <- elementIndex + 1
          lastIndex              <- k2+1
        }
      }
    }
  }
  elements[elementIndex] <- strtrimM(substr(input,lastIndex,nchar(input)))
  return(elements)
}


###############################################################
# strlocateall: Returns indices of occurrence of strings in strings
###############################################################
#' Returns indices of occurrence of strings in strings
#'
#' @param input A string to search in
#' @param searchString A string to search for
#' @return A list with the following components:
#' \item{start}{A vector with start indices}
#' \item{end}{A vector with end indices of matching}
#' @examples
#' strlocateall("abcdefgabcdefg","a")
#' strlocateall("abcdefgabcdefg","ab")
#' strlocateall("abcdefgabcdefg","x")
#' @export

strlocateall <- function(input,searchString) {
  x <- gregexpr(pattern=searchString,input,fixed=TRUE)
  start <- unlist(x)
  length <- attr(x[[1]],"match.length")
  end <- start+length-1
  if (length(start)==1 && start==-1)
    return(list(start=NULL,end=NULL))
  else
    return(list(start=start,end=end))
}

