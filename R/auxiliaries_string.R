###############################################################################
###############################################################################
# This file contains functions related to string manipulations
###############################################################################
###############################################################################


###############################################################################
# deString: converts a vector to a numeric vector
###############################################################################
# Converts vector to numeric vector
#
# A potentially non-numeric vector will be converted to a numeric vector.
# Numbers represented as strings will be converted to numerical values,
# non-numeric content will be conveted to NA.
#
# @param input A vector
# @return Numeric vector (non-numeric input elements converted to NA)
# @examples
# deString(c(1,2,3))
# deString(c('1',2,3))
# deString(c('1','a',3))

deString <- function(input) {
  suppressWarnings(as.numeric(input))
}


###############################################################
# strtrim: truncates strings by removing leading and trailing
#          white spaces (shadows the base strtrim function but
#          is more in line with strtrim in PHP and MATLAB)
###############################################################
# Truncates strings by removing leading and trailing white spaces
#
# Whitespace characters are the following: intToUtf8(c(9, 10, 11, 12, 13, 32))
#
# @param input A character string
# @return Character string with leading and trailing whitespaces removed
# @examples
# strtrim("abcde")
# strtrim("a b c d e")
# strtrim(" a b c d e    ")
# strtrim(intToUtf8(c(9, 10, 11, 12, 13, 32, 97,  32,  98,  32,  99,  32, 100,  32, 101, 32, 13)))

strtrim <- function(input) {
  outputText <- gsub("^\\s+|\\s+$", "", input)
  return(outputText)
}


###############################################################
# strrep: string replace function
###############################################################
# String replacement function
#
# @param origstr A character string in which to replace sub strings
# @param oldsubstr Search string to replace
# @param newsubstr A string to replace the search string if found
# @return Resulting output string
# @examples
# strrep("Hello World!","World","Kitty")
# strrep("a b c d e"," ","")

strrep <- function(origstr,oldsubstr,newsubstr) {
  outputText <- gsub(oldsubstr, newsubstr, origstr, fixed="TRUE")
  return(outputText)
}


###############################################################
# strremWhite: Removes all white spaces in strings
###############################################################
# Removes all white spaces in a string
#
# Whitespace characters are the following: intToUtf8(c(9, 10, 11, 12, 13, 32))
#
# @param input A character string
# @return Character string with whitespaces removed
# @examples
# strremWhite("abcde")
# strremWhite("a b c d e")
# strremWhite(" a b c d e    ")
# strremWhite(intToUtf8(c(9, 10, 11, 12, 13, 32, 97,  32,  98,  32,  99,  32, 100,  32, 101, 32, 13)))

strremWhite <- function(inputText) {
  outputText <- gsub("\\s+|\\s+", "", inputText)
  return(outputText)
}


###############################################################
# strmatch: Returns indices of elements in a vector of that match a search string
###############################################################
# Returns indices of elements in a vector of that match a search string
#
# Exact, case sensitive matching is used.
# Function can also be used for other objects than strings.
#
# @param inputVec A vector of strings
# @param searchString A string to search for (case sensitive,
#   exact matching) in inputVec
# @return A vector with indices of elements that have an exact match.
#   NULL is returned if string not found
# @examples
# strmatch(c("hello","test","a","e","hello"),"hello")
# strmatch(c("hello","test","a","e","hello"),"Hello")
# strmatch(c("hello","test","a","e","hello"),"xyz")
# strmatch(c(1,2,3,4,8,3),3)

strmatch <- function(inputVec,searchString) {
  yesno <- as.numeric(inputVec==searchString)
  if (sum(yesno)==0) return(NULL)
  indices <- yesno*c(1:length(yesno))
  return(indices[indices!=0])
}


###############################################################
# strexplode: Splits strings by separator
###############################################################
# Splits a string based on a separator substring
#
# By default a "," is used as separator.
#
# @param input A string to split
# @param separator A substring to use as separator
# @return A character vector with the split elements
# @examples
# strexplode(c("Hello,,Test,,X,,Z"))
# strexplode(c("Hello,,Test,,X,,Z"),",")
# strexplode(c("Hello,,Test,,X,,Z"),",,")
# strexplode(c("Hello,,Test,,X,,Z"),"§")

strexplode <- function(input,separator=",") {
  return(unlist(strsplit(input,separator)))
}

###############################################################
# strexplodePC: Splits strings by separator, ignore separator in parentheses
###############################################################
# Splits strings by separator only if separators not in parentheses
#
# By default a "," is used as separator.
#
# @param input A string to split
# @param separator A substring to use as separator
# @param group "round", "square", or "curly" allowing grouping by parantheses
#   of defined type in which elements are not exploded
# @return A character vector with the split elements
# @examples
# strexplodePC("(a,b),(c,d)")
# strexplodePC("((a,b),(c,d))") # compare with previous - lowest level considered
# strexplodePC("(a,b),(c,d))")
# strexplodePC("(a;b);(c;d)",";")

strexplodePC <- function(input,separator=",",group="round") {
  if (group=="round") expr <- "\\([^\\)]+\\)"
  if (group=="square") expr <- "\\[[^\\]]+\\]"
  if (group=="curly") expr <- "\\{[^\\}]+\\}"
  if (group!="round" && group!="square" && group!="curly") stop("strexplodePC: wrong group definition")

  # First find text inside parentheses
  m      <- gregexpr(expr,input,perl=TRUE)
  a      <- regmatches(input,m)
  a      <- unique(unlist(a))

  if (length(a) == 0)
    return(strexplode(input,separator))

  # Replace text inside parentheses
  y <- input
  for (k in 1:length(a))
    y <- strrep(y,a[k],paste("GROUPPIECE",k,sep=""))
  # Re-replace elements
  z <- strexplode(y,separator)
  for (k in 1:length(z))
    for (k2 in 1:length(a))
      z[k] <- strrep(z[k],paste("GROUPPIECE",k2,sep=""),a[k2])
  return(z)
}


###############################################################
# strlocateall: Returns indices of occurrence of strings in strings
###############################################################
# Returns indices of occurrence of strings in strings
#
# @param input A string to search in
# @param searchString A string to search for
# @return A list with the following components:
# \item{start}{A vector with start indices}
# \item{end}{A vector with end indices of matching}
# @examples
# strlocateall("abcdefgabcdefg","a")
# strlocateall("abcdefgabcdefg","ab")
# strlocateall("abcdefgabcdefg","x")

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

