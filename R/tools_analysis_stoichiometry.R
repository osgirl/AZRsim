###############################################################################
###############################################################################
# This file contains functions for stoichiometric matrix determination
###############################################################################
###############################################################################

###############################################################################
# AZRstoichiometry: function returing the stoichiometric matrix for a an AZRmodel
###############################################################################
#' Determination of stoichiometric matrix
#'
#' Determines the stoichiometric matrix for the given AZRmodel.
#'
#' In order determine the stoichiometric matrix, the differential
#' equations for the components have to be
#' expressed in terms of reaction rates. The stoichiometric constants need
#' to be numeric and multiplied to the reaction terms.
#'
#' Example: d/dt(A) = -1*Re1+2*Re3-Re5+3.141*Re8
#'
#' Especially when importing models from SBML the right hand side of the
#' ODEs might show a correction term that is needed for transport between two
#' different compartments in the case that the species is defined in
#' concentration units. In this case the ODE can look as follows:
#'
#'    d/dt(B) = (-1*Re1+2*Re3-Re5+3.141*Re8)/compartmentsize
#'
#' This syntax is also accepted. In this case the stoichiometric elements
#' in the parenthesis will be divided by 'compartmentsize'.
#'
#' The 'compartmentsize' is only allowed to be a parameter! It can not be a
#' numeric value, a state, a variable, or a function!
#'
#' Except the above shown pair of parentheses, no additional parentheses are
#' allowed to appear in the ODE definitions. Otherwise, and in case that not
#' only reaction rates are present in the ODE expression, the corresponding
#' component is not taken into account for the stoichiometric matrix.
#'
#' @param model An AZRmodel
#' @param raw Flag used to force the AZRstoichiometry function
#'  not to correct the elements of the stoichiometric matrix for the
#'  compartment information in cases where species are given in
#'  concentrations. This is needed for model construction (BC type of
#'  representation and for other things). This flag should be kept on "FALSE"
#'  for the casual user of this function.
#' @return A list with the following components:
#' \item{N}{The stoichimoetric matrix}
#' \item{statenames}{Vector with names of states in N}
#' \item{reacnames}{Vector with reaction names in N}
#' \item{reacreversible}{Vector with reversible flags for the reactions in reacnames}
#' @examples
#' filename <- system.file(package="AZRsim","examples","NovakTyson.txt")
#' model <- AZRmodel(filename)
#' AZRstoichiometry(model)
#' @export

AZRstoichiometry <- function (model,raw=TRUE) {

  if (!is.AZRmodel(model))
    stop("AZRstoichiometry: model argument is not an AZRmodel")

  # Get state, parameter, reaction information
  stateInfo <- getAllStatesAZRmodel(model)
  paramInfo <- getAllParametersAZRmodel(model)
  reacInfo <- getAllReactionsAZRmodel(model)

  # Return in trivial case (no reactions present)
  if (getNumberOfReactionsAZRmodel(model)==0)
    return(list(N=NULL, statenames=NULL, reacnames=NULL, reacreversible=NULL))

  ###############################################
  # CHECK ODEs for reactions
  ###############################################
  # check if given ODE expression contains only reaction terms - then return
  # the stoichiometry for this ODE. Otherwise if not only reactions are
  # present, return an empty vector.
  # in cases that the reaction rate is adjusted by the compartment volume
  # (happens especially in cases of import of SBML models containing species
  # with concentration rates and compartment volumes different from one) the
  # adjustment term and the needed parentheses should be detected and
  # neglected.
  getStoichiometryInformation <- function(ODE) {
    ODE <- strremWhite(ODE)
    # check if a scaling by the compartment volume is done
    # in this case the expected syntax is
    # ODE = ("reactionterms")/compartmentvolume
    numberOpenParentheses <- length(strlocateall(ODE,"(")$start)
    numberClosedParentheses <- length(strlocateall(ODE,")")$start)
    compartmentSize <- 1
    if (numberOpenParentheses != numberClosedParentheses)
      stop("getStoichiometryInformation: parentheses not in pairs")
    if (numberOpenParentheses > 1)
      return(NULL) # only one pair of parentheses allowed
    if (numberOpenParentheses==1) {
      if (substr(ODE,1,1) != "(")
        # if a parenthesis exists in the ODE it has to be the first character
        return(NULL)

      # One parenthesis present and it is the first char in the ODE
      # assume that this is due to a adjustement to compartment sizes
      # cut out the content of the parentheses
      closePar <- strlocateall(ODE,")")
      ODEinpar <- substr(ODE,2,closePar$start-1)
      rest <- substr(ODE,closePar$start+1,nchar(ODE))
      # first character needs to be a '/' based on assumption of compartment scaling
      if (substr(rest,1,1) != "/")
        # if not then return
        return(NULL)

      rest <- substr(rest,2,nchar(rest))

      # check if this rest corresponds to a parameter name and if yes get its value
      index <- strmatch(paramInfo$paramnames,rest)
      if (is.null(index))
        # not a parameter name => return
        return(NULL)

      # get compartment size as value of parameter
      compartmentSize = paramInfo$paramvalues[index]

      ODE <- ODEinpar
    }

    ODE <- paste(ODE,"+",sep="")
    # first explode the ODE in additive terms
    terms <- c()
    termIndex <- 1
    # check the sign of the first term (first character in string)
    if (substr(ODE,1,1)=="-") {
      signCurrent <- -1
      lastIndex <- 2
    } else {
      if (substr(ODE,1,1)=="+") {
        signCurrent <- +1
        lastIndex <- 2
      } else {
        signCurrent = +1
        lastIndex = 1
      }
    }

    # explode in terms, check if term has the right format and if the
    # second term features a reactionname. then construct the row of the
    # stoichiometric matrix
    Nrow <- as.vector(matrix(0,1,length(reacInfo$reacnames)))
    startk <- lastIndex
    for (k in startk:nchar(ODE)) {
      # check for positive or negative term
      if (substr(ODE,k,k) == '+' || substr(ODE,k,k) == '-') {
        element <- substr(ODE,lastIndex,k-1)
        # check the element if composed of term in the right format
        multIndex <- strlocateall(element,"*")$start
        if (is.null(multIndex)) {
          stoichiometry <- signCurrent*1
          reactionterm <- element
        } else {
          if (length(multIndex) > 1)
            # to many multiplication signs (only one allowed)
            return(NULL)

          # only one multiplication and factor needs to be numeric (otherwise
          # absStoichiometry is NA and later return with NULL)
          suppressWarnings(absStoichiometry <- as.numeric(substr(element,1,multIndex-1)))
          stoichiometry <- signCurrent*absStoichiometry
          reactionterm = substr(element,multIndex+1,nchar(element))
        }

        # find the index of the reaction name and add the
        # stoichiometric information to Nrow
        indexReaction <- strmatch(reacInfo$reacnames,reactionterm)
        if (is.null(indexReaction))
          return(NULL)
        if(is.na(stoichiometry))
          return(NULL)

        Nrow[indexReaction] <- Nrow[indexReaction] + stoichiometry

        # increment
        termIndex <- termIndex + 1
        lastIndex <- k+1
        if (substr(ODE,k,k) == "+")
          signCurrent <- +1
        else
          signCurrent <- -1
      }
    }

    # adjust stoichiometries with compartment size
    if (!raw)
      Nrow <- Nrow / compartmentSize

    return(Nrow)
  }

  ###########################################################
  # Go through all ODEs and check for which states the ODEs only
  # consist of reaction terms - this returns already positive and negative
  # terms and coefficients
  ###########################################################
  N <- c()
  statenames <- c()

  for (k in 1:getNumberOfStatesAZRmodel(model)) {
    # check if ODE contains only reaction terms - otherwise do not
    # consider the current state as component for the stoichiometric matrix
    Nrow <- getStoichiometryInformation(model$states[[k]]$ODE)
    if (!is.null(Nrow) && !(NA %in% Nrow)) {
      N <- rbind(N,Nrow)
      statenames = c(statenames,stateInfo$statename[k])
    }
  }

  return(list(N=N,statenames=statenames,reacnames=reacInfo$reacnames,reacreversible=reacInfo$reacreversible))
}
