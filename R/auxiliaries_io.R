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
fileparts <- function(filename.with.path){
  pathname <- dirname(filename.with.path)
  filename <- basename(filename.with.path)
  fileext <- gsub(".*(\\.[^\\.]*)$","\\1",filename)
  filename <- gsub("(.*)(\\.[^\\.]*)$","\\1",filename)
  if(fileext==filename) fileext <- ""
  return(list(pathname=pathname,filename=filename,fileext=fileext))
}


#' Writes a formatted string to a text file
#'
#' @details
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
filewrite <- function(text,filename) {
  fid <- fopen(filename, mode="w")
  write(text, fid)
  fclose(fid)
}


#' Reads a text file and returns it as string
#'
#' @details
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
fileread <- function(filename,collapserows=TRUE) {
  fid <- fopen(filename, mode="r")
  text <- readLines(fid)
  fclose(fid)
  if (collapserows) {
    text <- paste(text,collapse="\n")
  }
  return(text)
}

#' Creates a folder if it does not yet exist
#'
#' @param pathdir String with folder path
#' @return None
#' @examples
#' \dontrun{
#' mkdirp("test/test")
#' mkdirp("../test")
#' }
#' @export
mkdirp <- function(pathdir) {
  if (!file.exists(pathdir)) dir.create(pathdir,recursive='TRUE')
}


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
fopen <- function(filename,mode="w") {
  if (mode=="w") mkdir(fileparts(filename)$pathname)
  fid <- file(filename,open=mode)
  return(fid)
}


#' Closes a file
#'
#' @param fid File ID
#' @return None
#' @examples
#' \dontrun{
#' fclose(fid)
#' }
fclose <- function(fid) {
  close(fid)
}


#' Writes text into an opened file
#'
#' @details
#' File needs to be opened with fopen and closed with fclose
#'
#' @param fid File ID
#' @param text Text to write - can be formatted with escaped chars
#' @return None
#' @examples
#' \dontrun{
#' fwrite(fid,text)
#' }
fwrite <- function(fid,text) {
  write(text,fid)
}
