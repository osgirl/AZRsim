###############################################################################
###############################################################################
###############################################################################
###############################################################################

###############################################################################
# .onAttach
###############################################################################
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("AZRsim: Powerful dynamic model representation and simulation",
                              "        Part of the 'AZRtools Suite of R packages'",
                              "        Enjoy!",sep="\n"))
}

###############################################################################
# .onLoad
###############################################################################
.onLoad <- function(libname, pkgname) {
}

###############################################################################
# .onUnload
###############################################################################
.onUnload <- function(libpath) {
}

###############################################################################
# Halleluja CRAN check ...
# globalVariables is a hideous hack and I will never use it. – hadley Sep 25 '12 at 16:10
# @hadley you shouldn't say you'll never use things when only two years later you think it's fine – hadley Dec 31 '14 at 19:31
# http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
###############################################################################
globalVariables(c("TIME", "INPUT", "TIME_DOSE_EFFECT_START", "time"))
