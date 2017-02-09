###############################################################################
###############################################################################
# This file contains functions related to (file) input/output
###############################################################################
###############################################################################

###############################################################################
# fileparts: writes a formatted string to a text file
###############################################################################
#' fileparts function
#'
#' Get the various parts of a file with path string.
#'
#' @param filename.with.path A string of a filename with a path
#'
#' @return A list with the following components:
#' \item{pathname}{The path name}
#' \item{filename}{The file name}
#' \item{fileext}{The file extension}
#'
#' @examples
#' fileparts("a/b/c/d.R")
#' fileparts("e")
#' fileparts("f.R")
#' @export

fileparts <- function(filename.with.path){
  pathname <- dirname(filename.with.path)
  filename <- basename(filename.with.path)
  fileext <- gsub(".*(\\.[^\\.]*)$","\\1",filename)
  filename <- gsub("(.*)(\\.[^\\.]*)$","\\1",filename)
  if(fileext==filename) fileext <- ""
  return(list(pathname=pathname,filename=filename,fileext=fileext))
}


###############################################################################
# filewrite: writes a formatted string to a text file
###############################################################################
#' Writes a formatted string to a text file
#'
#' This is a wrapper to the write function. Ideally, use it only for text to
#' write to a file. If the file with the filename already exists, it will be
#' overwritten. If a path is provided and the folder does not exist, the folder
#' will be created first.
#'
#' @param text A character string
#' @param filename, possibly including the path
#' @return None
#' @examples
#' \dontrun{
#' filewrite("Hello World!","filename.txt")
#' filewrite("Hello World!","test/filename.txt")
#' }
#' @export

filewrite <- function(text,filename) {
  fid <- fopen(filename, mode="w")
  write(text, fid)
  fclose(fid)
}


###############################################################################
# fileread: reads contents of a text file
###############################################################################
#' Reads a text file and returns it as string
#'
#' Should only be used on text files. Will read line by line until the end and
#' return the results in a character string, using "\\n" as separator between
#' read lines.
#'
#' @param filename The filename, possibly including the path
#' @param collapserows FALSE: each row of the txt file is returned as a separate list element
#'                     TRUE: a single block is returned for the whole text
#' @return Character string with files contents
#' @examples
#' \dontrun{
#' filewrite("Hello World!","filename.txt")
#' fileread("filename.txt")
#' filewrite("Hello\nWorld!","filename.txt")
#' fileread("filename.txt")
#' fileread("filename.txt",collapserows=FALSE)
#' }
#' @export

fileread <- function(filename,collapserows=TRUE) {
  fid <- fopen(filename, mode="r")
  text <- readLines(fid)
  fclose(fid)
  if (collapserows) {
    text <- paste(text,collapse="\n")
  }
  return(text)
}

###############################################################################
# mkdir: Creates a folder if it does not yet exist
###############################################################################
#' Creates a folder if it does not yet exist
#'
#' @param pathdir String with folder path
#' @return None
#' @examples
#' \dontrun{
#' mkdir("test/test")
#' mkdir("../test")
#' }
#' @export

mkdir <- function(pathdir) {
  if (!file.exists(pathdir)) dir.create(pathdir,recursive='TRUE')
}


###############################################################################
# rmdir: Removes a folder
###############################################################################
#' Removes a folder
#'
#' @param pathdir String with folder path
#' @return None
#' @examples
#' \dontrun{
#' rmdir("test/test")
#' rmdir("../test")
#' }
#' @export

rmdir <- function(pathdir) {
  unlink(pathdir,recursive = 'TRUE')
}


###############################################################################
# fopen: Opens a file and creates the folder if needed
###############################################################################
#' Opens a file and creates the folder if needed
#'
#' @param filename String with file, possibly including path
#' @param mode "w" for writing, "r" for reading
#' @return File ID
#' @examples
#' \dontrun{
#' fopen("e:/test/test.txt")
#' fopen("../test.R")
#' }
#' @export

fopen <- function(filename,mode="w") {
  if (mode=="w") mkdir(fileparts(filename)$pathname)
  fid <- file(filename,open=mode)
  return(fid)
}


###############################################################################
# fclose: Closes a file
###############################################################################
#' Closes a file
#'
#' @param fid File ID
#' @return None
#' @examples
#' \dontrun{
#' fclose(fid)
#' }
#' @export

fclose <- function(fid) {
  close(fid)
}


###############################################################################
# fwrite: Writes text into an opened file
###############################################################################
#' Writes text into an opened file
#'
#' File needs to be opened with fopen and closed with fclose
#'
#' @param fid File ID
#' @param text Text to write - can be formatted with escaped chars
#' @return None
#' @examples
#' \dontrun{
#' fwrite(fid,text)
#' }
#' @export

fwrite <- function(fid,text) {
  write(text,fid)
}
