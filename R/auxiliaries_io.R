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
