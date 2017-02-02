###############################################################################
###############################################################################
# This file contains the AZRmodel function and simple methods for handling AZRmodels
###############################################################################
###############################################################################


###############################################################################
# AZRmodel: function returing S3 object "AZRmodel"
###############################################################################
#' AZRmodel class
#'
#' AZRmodel objects represent dynamic, ODE based models that can be simulated.
#' The class constructor function AZRmodel() can create such objects from
#' text file descriptions of these models. If no input argument is provided,
#' an empty AZRmodel object is returned.
#'
#' @param input A string with the full filename (including absolute or
#'   relative path) to a text file, describing an AZRmodel.
#' @param simFlag If set to TRUE (default) then the model will be prepared for
#'   simulations (dosing, constraints, etc.) and all required sim functions will
#'   be generated and added to the AZRmodel attributes. Choose this option when
#'   you want to simulate. In cases were the model structure and the information
#'   about inputs is important, then choose "FALSE".
#' @return An AZRmodel object. If input is NULL, an empty object is returned.
#' @examples
#' AZRmodel()
#' filename <- system.file(package="AZRsim","examples","NovakTyson.txt")
#' AZRmodel(filename)
#' @export

AZRmodel <- function (input=NULL,simFlag=TRUE) {

  #################################
  # Handle optional input arguments
  #################################
  # No optional input arguments present

  #################################
  # Checking of input arguments
  #################################
  if(!is.null(input) && !file.exists(input)) {
    # input does not point to an existing file on the file system => error
    stop("AZRmodel: Provided input argument does not point to a file on the filesystem.")
  }

  #################################
  # Initialize empty AZRmodel
  #################################
  model <- createEmptyAZRmodel()

  #################################
  # Return since no input model file defined
  #################################
  if (is.null(input)) {
    return(model)
  }

  #################################
  # Get extension of model file
  #################################
  fileInfo <- fileparts(input)

  #################################
  # Try to import the file
  #################################
  if (fileInfo$fileext==".txt") {
    model <- importTxtAZRmodel(model,input)
  } else {
    if (fileInfo$fileext==".txtbc") {
      model <- importTxtBcAZRmodel(model,input)
    } else {
      if (fileInfo$fileext==".xml") {
        stop("Import of AZRmodels from SBML not supported in AZR Tools!")
      } else {
        stop("Unknown model file extension.")
      }
    }
  }

  #################################
  # Check names of components
  #################################
  checkNamesAZRmodel(model)


  #################################
  # Add the original model as attribute (before generation of simulation functions
  # and required update of the model
  #################################
  attr(model,"originalModel") <- model

  #################################
  # Update the models math to for simulation purposes
  # Constraints
  # Inputs / Dosing
  # Namespace for deSolve (only in simfunctions)
  #################################
  if (simFlag) {
    model <- genSimFunctions(model)
  }

  #################################
  # Return model
  #################################
  # construct the model object
  class(model) <- "AZRmodel"
  return(model)
}


###############################################################################
# AZRexportAZRmodel: exports an AZRmodel as .txt or .txtbc file
###############################################################################
#' Export of AZRmodel
#'
#' Export of an AZRmodel to either an ODE based .txt or a biochemical reaction
#' based .txtbc file. Note that .txtbc files also might contain ODEs, as not all
#' models allow for the biochemical reaction notation.
#'
#' @param model An AZRmodel to be exported.
#' @param filename Full path with filename to export the model
#' @param useBC Flag for the use of the biochemical notation if TRUE
#' @param useSIM Flag to export the model with potential changes for simulation
#'   purposes if TRUE. If FALSE, the originally imported model will be exported.
#' @return None
#' @examples
#' model <- AZRmodel()
#' AZRexportAZRmodel(model)
#' filename <- system.file(package="AZRsim","examples","NovakTyson.txt")
#' model <- AZRmodel(filename)
#' AZRexportAZRmodel(model,"filename")
#' AZRexportAZRmodel(model,'filename',useBC=TRUE)
#' @export

AZRexportAZRmodel <- function (model, filename=NULL, useBC=FALSE, useSIM=FALSE) {

  if (!is.AZRmodel(model))
    stop("AZRexportAZRmodel: input argument is not an AZRmodel")

  if (!useSIM && !is.null(attr(model,"originalModel"))) model <- attr(model,"originalModel")

  # Convert to either TEXT or TEXTBC
  if (!useBC)
    exportTxtAZRmodel(model,filename)
  else
    exportTxtBcAZRmodel(model,filename)
}


###############################################################################
# createEmptyAZRmodel
###############################################################################
# Creates an empty AZRmodel object
#
# @return An empty AZRmodel object.
# @examples
# createEmptyAZRmodel()

createEmptyAZRmodel <- function () {

  model                     <- list()
  model$name                <- "unnamed_model"
  model$notes               <- "model notes"
  model$states              <- list()
  model$algebraic           <- list()
  model$parameters          <- list()
  model$variables           <- list()
  model$reactions           <- list()
  model$events              <- list()
  model$functions           <- list()
  model$inputs              <- list()
  model$outputs             <- list()

  # construct the model object
  class(model) <- "AZRmodel"
  return(model)
}


###############################################################################
# Generic function overload (is) - need to be exported
###############################################################################
#' Check function if it is an AZRmodel
#'
#' @param input AZRmodel object
#' @return TRUE or FALSE
#' @examples
#'   is.AZRmodel(AZRmodel())
#' @export
is.AZRmodel <- function(input) {
  methods::is(input,"AZRmodel")
}

###############################################################################
# Generic function overload (print) - need to be exported
###############################################################################
#' Generic print function for AZRmodels
#'
#' @param x AZRmodel object
#' @param ... Additional unused arguments
#' @examples
#'   print(AZRmodel())
#' @export
print.AZRmodel <- function(x, ...) {
  cat("\tAZRmodel\n\t========\n")
  cat("\tName:                      ", x$name,"\n")
  cat("\tNumber States:             ", getNumberOfStatesAZRmodel(x),"\n")
  if (getNumberOfAlgebraicAZRmodel(x)>0)
    cat("\tNumber Algebraic States:   ", getNumberOfAlgebraicAZRmodel(x),"\n")
  cat("\tNumber Parameters:         ", getNumberOfParametersAZRmodel(x),"\n")
  cat("\tNumber Variables:          ", getNumberOfVariablesAZRmodel(x),"\n")
  cat("\tNumber Reactions:          ", getNumberOfReactionsAZRmodel(x),"\n")
  cat("\tNumber Functions:          ", getNumberOfFunctionsAZRmodel(x),"\n")
  if (getNumberOfEventsAZRmodel(x)>0)
    cat("\tNumber Events:             ", getNumberOfEventsAZRmodel(x),"\n")
  if (getNumberOfInputsAZRmodel(x)>0)
    cat("\tNumber Inputs:             ", getNumberOfInputsAZRmodel(x),"\n")
  if (getNumberOfOutputsAZRmodel(x)>0)
    cat("\tNumber Outputs:            ", getNumberOfOutputsAZRmodel(x),"\n")
  if (getNumberOfAlgebraicAZRmodel(x) > 0)
    cat("\tAlgebraic definitions of states present in the model. Such states
        can not be handled by other functions. Please check if these states
        are really required. Many modeling tools (e.g. simbiology) add these
        for simple moiety conservations, that are better represented as
        variables. If you still need them, contact us (info@intiquan.com
        and we will see what we can do! In the meantime, consider the use
        of IQM Tools (http://www.intiquan.com/iqm-tools/)\n")
  if (!hasonlynumericICsAZRmodel(x))
    cat("\tNon-numeric initial conditions are present in the model.\n")
  if (hasfastreactionsAZRmodel(x))
    cat("\tFast reactions defined in the model. Such reactions can not be
        handled by other functions. Please check if these reactions
        are really required. If you need them, contact us (info@intiquan.com
        and we will see what we can do! In the meantime, consider the use
        of IQM Tools (http://www.intiquan.com/iqm-tools/)\n")
}


###############################################################################
# checkNamesAZRmodel
###############################################################################
# Checks names of model components and issues warnings or errors
checkNamesAZRmodel <- function(model) {
  stateNames <- toupper(getAllStatesAZRmodel(model)$statenames)
  paramNames <- toupper(getAllParametersAZRmodel(model)$paramnames)
  varNames   <- toupper(getAllVariablesAZRmodel(model)$varnames)
  reacNames  <- toupper(getAllReactionsAZRmodel(model)$reacnames)

  # Check for single char name (should be avoided)
  if (sum(as.numeric(nchar(stateNames) == 1)) > 0)
    warning("checkNamesAZRmodel: AZRmodel contains state names with a single character name. Try to avoid that if you plan to use NONMEM or MONOLIX")
  if (sum(as.numeric(nchar(paramNames) == 1)) > 0)
    warning("checkNamesAZRmodel: AZRmodel contains parameter names with a single character name. Try to avoid that if you plan to use NONMEM or MONOLIX")
  if (sum(as.numeric(nchar(varNames) == 1)) > 0)
    warning("checkNamesAZRmodel: AZRmodel contains variable names with a single character name. Try to avoid that if you plan to use NONMEM or MONOLIX")
  if (sum(as.numeric(nchar(reacNames) == 1)) > 0)
    warning("checkNamesAZRmodel: AZRmodel contains reaction names with a single character name. Try to avoid that if you plan to use NONMEM or MONOLIX")

  # Define reserved words - case insensitive
  reservedWords <- toupper(c("F","G","H","gt","ge","lt","le","mod","and","or","multiply",
                             "piecewise","interp0","interp1","interpcs","interpcse"))

  # Check if model contains elements with names matching reserved words
  if (length(intersect(stateNames,reservedWords)) > 0)
    stop(paste("checkNamesAZRmodel: model contains the following state name(s) that is(are) reserved word(s): ",intersect(stateNames,reservedWords),sep=""))
  if (length(intersect(paramNames,reservedWords)) > 0)
    stop(paste("checkNamesAZRmodel: model contains the following parameter name(s) that is(are) reserved word(s): ",intersect(paramNames,reservedWords),sep=""))
  if (length(intersect(varNames,reservedWords)) > 0)
    stop(paste("checkNamesAZRmodel: model contains the following variable name(s) that is(are) reserved word(s): ",intersect(varNames,reservedWords),sep=""))
  if (length(intersect(reacNames,reservedWords)) > 0)
    stop(paste("checkNamesAZRmodel: model contains the following reaction name(s) that is(are) reserved word(s): ",intersect(reacNames,reservedWords),sep=""))
}


###############################################################################
# hasonlynumericICsAZRmodel
###############################################################################
# Checks if the model contains only numeric initial conditions
#
# @param model An AZRmodel object
# @return TRUE if model contains non-numerical initial conditions, FALSE otherwise
# @examples
# model <- exNovakTyson
# hasonlynumericICsAZRmodel(model)
hasonlynumericICsAZRmodel <- function(model) {

  if (!is.AZRmodel(model))
    stop("hasonlynumericICsAZRmodel: input argument is not an AZRmodel")

  if (getNumberOfStatesAZRmodel(model) < 1) return(TRUE)

  for (k in 1:getNumberOfStatesAZRmodel(model))
    if (!isnumericVector(model$state[[k]]$IC))
      return(FALSE)

  return(TRUE)
}

###############################################################################
# hasfastreactionsAZRmodel
###############################################################################
# Checks if the model contains fast reactions
#
# @param model An AZRmodel object
# @return TRUE if model contains fast reactions, FALSE otherwise
# @examples
# model <- exNovakTyson
# hasfastreactionsAZRmodel(model)
hasfastreactionsAZRmodel <- function(model) {

  if (!is.AZRmodel(model))
    stop("hasfastreactionsAZRmodel: input argument is not an AZRmodel")

  if (getNumberOfReactionsAZRmodel(model) < 1) return(FALSE)

  for (k in 1:getNumberOfReactionsAZRmodel(model))
    if (model$reactions[[k]]$fast)
      return(TRUE)

  return(FALSE)
}

###############################################################################
# hasalgebraicAZRmodel
###############################################################################
# Checks if the model contains algebraic states
#
# @param model An AZRmodel object
# @return TRUE if model contains algebraic states, FALSE otherwise
# @examples
# model <- exNovakTyson
# hasalgebraicAZRmodel(model)
hasalgebraicAZRmodel <- function(model) {

  if (!is.AZRmodel(model))
    stop("hasalgebraicAZRmodel: input argument is not an AZRmodel")

  if (getNumberOfAlgebraicAZRmodel(model) < 1) return(FALSE)

  return(TRUE)
}

###############################################################################
# Functions to get numbers of model elements
###############################################################################

# Returns number of states in an AZRmodel
#
# @param model An AZRmodel object
# @return Number of states
# @export
getNumberOfStatesAZRmodel <- function(model) {
  length(model$states)
}

# Returns number of algebraic states in an AZRmodel
#
# @param model An AZRmodel object
# @return Number of algebraic states
# @export
getNumberOfAlgebraicAZRmodel <- function(model) {
  length(model$algebraic)
}

# Returns number of parameters in an AZRmodel
#
# @param model An AZRmodel object
# @return Number of parameters
# @export
getNumberOfParametersAZRmodel <- function(model) {
  length(model$parameters)
}

# Returns number of variables in an AZRmodel
#
# @param model An AZRmodel object
# @return Number of variables
# @export
getNumberOfVariablesAZRmodel <- function(model) {
  length(model$variables)
}

# Returns number of reactions in an AZRmodel
#
# @param model An AZRmodel object
# @return Number of reactions
# @export
getNumberOfReactionsAZRmodel <- function(model) {
  length(model$reactions)
}

# Returns number of functions in an AZRmodel
#
# @param model An AZRmodel object
# @return Number of functions
# @export
getNumberOfFunctionsAZRmodel <- function(model) {
  length(model$functions)
}

# Returns number of events in an AZRmodel
#
# @param model An AZRmodel object
# @return Number of events
# @export
getNumberOfEventsAZRmodel <- function(model) {
  length(model$events)
}

# Returns number of inputs in an AZRmodel
#
# @param model An AZRmodel object
# @return Number of inputs
# @export
getNumberOfInputsAZRmodel <- function(model) {
  length(model$inputs)
}

# Returns number of outputs in an AZRmodel
#
# @param model An AZRmodel object
# @return Number of outputs
# @export
getNumberOfOutputsAZRmodel <- function(model) {
  length(model$outputs)
}

# Returns number of Event Assignments in an Event in an AZRmodel
#
# @param model An AZRmodel object
# @param eventindex Index of the event to check
# @return Number of event assignments
# @export
getNumberOfEventassignmentsAZRmodel <- function(model,eventindex) {
  if (eventindex > getNumberOfEventsAZRmodel(model))
    stop("getNumberOfEventassignmentsAZRmodel: eventindex larger than the number of events in the model.")
 length(model$events[[eventindex]]$assignment)
}


###############################################################################
# State handling functions
# addStateAZRmodel
# getStateAZRmodel
# setStateAZRmodel
# delStateAZRmodel
###############################################################################

# Add a new state to an AZRmodel
#
# @param model An AZRmodel
# @param name String with state name
# @param IC Initial condition (numeric or string)
# @param ODE String with RHS of ODE
# @param lowConstraint Value for lower bound of state variable
# @param highConstraint Value for upper bound of state variable
# @param type String with type of SBML element (isSpecie, isParameter, isCompartment)
# @param compartment String with name of SBML compartment for species or outside compartment for compartment)
# @param unittype String with type of unit (concentration, amount)
# @param notes String with notes about the state
# @return An AZRmodel object with appended state variable
# @examples
# model <- AZRmodel()
# model <- addStateAZRmodel(model,'Cyclin',IC=0.3,'Re1-Re2')
# model <- addStateAZRmodel(model,'Parasites',IC=1e9,'(GR-KR)*Paasites')
# @export
addStateAZRmodel <- function(model,
                             name = NULL,
                             IC = NULL,
                             ODE = NULL,
                             lowConstraint = NULL,
                             highConstraint = NULL,
                             type = NULL,
                             compartment = NULL,
                             unittype = NULL,
                             notes = NULL) {

  if (!is.AZRmodel(model))
    stop("addStateAZRmodel: model argument is not an AZRmodel")

  if (is.null(name))
    stop("addStateAZRmodel: name is a required input argument")

  if (is.null(IC))
    stop("addStateAZRmodel: IC is a required input argument")

  if (is.null(ODE))
    stop("addStateAZRmodel: ODE is a required input argument")

  if (is.null(lowConstraint) && !is.null(highConstraint))
    stop("addStateAZRmodel: If highConstraint is defined also lowConstraint needs to be defined")

  if (!is.null(lowConstraint) && is.null(highConstraint))
    stop("addStateAZRmodel: If lowConstraint is defined also highConstraint needs to be defined")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("addStateAZRmodel: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("addStateAZRmodel: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

  stateInfo = list(name = name, IC = IC, ODE = ODE,
                   lowConstraint = lowConstraint,
                   highConstraint = highConstraint,
                   type = type,
                   compartment = compartment,
                   unittype = unittype,
                   notes = notes)

  model$states[[length(model$states)+1]] <- stateInfo
  return(model)
}


# Get AZRmodel state information
#
# @param model An AZRmodel
# @param index Index of the state in the model
# @return A list with the state information
# @examples
# model <- exNovakTyson
# getStateAZRmodel(model,1)
# @export
getStateAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("getStateAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfStatesAZRmodel(model) < index)
    stop("getStateAZRmodel: value of index larger than number of states.")

  if (index < 1)
    stop("getStateAZRmodel: value of index should be larger than 0.")

  return(model$states[[index]])
}


# Set state information in an AZRmodel
#
# Only provided information is updated in the state definition.
#
# @param model An AZRmodel
# @param index Numerical index of the state to update
# @param name String with state name
# @param IC Initial condition (numeric or string)
# @param ODE String with RHS of ODE
# @param lowConstraint Value for lower bound of state variable
# @param highConstraint Value for upper bound of state variable
# @param type String with type of SBML element (isSpecie, isParameter, isCompartment)
# @param compartment String with name of SBML compartment for species or outside compartment for compartment)
# @param unittype String with type of unit of SBML species (concentration, amount)
# @param notes String with notes about the state
# @return An AZRmodel object with updated state variable
# @examples
# model <- exNovakTyson
# model <- setStateAZRmodel(model,1,IC=0.5)
# @export
setStateAZRmodel <- function(model, index,
                             name = NULL,
                             IC = NULL,
                             ODE = NULL,
                             lowConstraint = NULL,
                             highConstraint = NULL,
                             type = NULL,
                             compartment = NULL,
                             unittype = NULL,
                             notes = NULL) {

  if (!is.AZRmodel(model))
    stop("setStateAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfStatesAZRmodel(model) < index)
    stop("setStateAZRmodel: value of index larger than number of states.")

  if (index < 1)
    stop("setStateAZRmodel: value of index should be larger than 0.")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("setStateAZRmodel: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("setStateAZRmodel: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

  if (!is.null(name))           model$states[[index]]$name           <- name
  if (!is.null(IC))             model$states[[index]]$IC             <- IC
  if (!is.null(ODE))            model$states[[index]]$ODE            <- ODE
  if (!is.null(lowConstraint))  model$states[[index]]$lowConstraint  <- lowConstraint
  if (!is.null(highConstraint)) model$states[[index]]$highConstraint <- highConstraint
  if (!is.null(type))           model$states[[index]]$type           <- type
  if (!is.null(compartment))    model$states[[index]]$compartment    <- compartment
  if (!is.null(unittype))       model$states[[index]]$unittype       <- unittype
  if (!is.null(notes))          model$states[[index]]$notes          <- notes

  return(model)
}


# Delete AZRmodel state
#
# @param model An AZRmodel
# @param index Index of the state in the model
# @return An AZRmodel with the indexed state removed
# @examples
# model <- exNovakTyson
# delStateAZRmodel(model,1)
# @export
delStateAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("delStateAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfStatesAZRmodel(model) < index)
    stop("delStateAZRmodel: value of index larger than number of states.")

  if (index < 1)
    stop("delStateAZRmodel: value of index should be larger than 0.")

  # check if state has inputs
  if (getNumberOfInputsAZRmodel(model) > 0) {
    inputstates = c()
    for (k in 1:getNumberOfInputsAZRmodel(model)) inputstates = c(inputstates,model$inputs[[k]]$stateindex)
    if (index %in% inputstates)
      stop("delStateAZRmodel: The state contains inputs. Please delete them first and then the state.")
  }

  model$states[[index]] <- NULL

  return(model)
}


###############################################################################
# Parameter handling functions
# addParameterAZRmodel
# getParameterAZRmodel
# setParameterAZRmodel
# delParameterAZRmodel
###############################################################################

# Add a new parameter to an AZRmodel
#
# @param model An AZRmodel
# @param name String with parameter name
# @param value Numerical value of parameter
# @param type String with type of SBML element (isSpecie, isParameter, isCompartment)
# @param compartment String with name of SBML compartment for species or outside compartment for compartment)
# @param unittype String with type of unit of SBML species (concentration, amount)
# @param notes String with notes about the parameter
# @param estimate Flag for parameter being estimated
# @param regressor Flag for parameter being a regressor in estimation
# @return An AZRmodel object with appended parameter
# @examples
# model <- AZRmodel()
# addParameterAZRmodel(model,'ka',2,notes="Hello World")
# addParameterAZRmodel(model,'F',value=0.5)
# @export
addParameterAZRmodel <- function(model,
                             name = NULL,
                             value = NULL,
                             type = NULL,
                             compartment = NULL,
                             unittype = NULL,
                             notes = NULL,
                             estimate = FALSE,
                             regressor = FALSE) {

  if (!is.AZRmodel(model))
    stop("addParameterAZRmodel: model argument is not an AZRmodel")

  if (is.null(name))
    stop("addParameterAZRmodel: name is a required input argument")

  if (is.null(value))
    stop("addParameterAZRmodel: value is a required input argument")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("addParameterAZRmodel: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("addParameterAZRmodel: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

  if (estimate && regressor)
    stop("addStateAZRmodel: estimate and regressor can not be TRUE at the same time")

  paramInfo <- list(name=name, value=value,
                   type = type,
                   compartment = compartment,
                   unittype = unittype,
                   notes = notes,
                   estimate = estimate,
                   regressor = regressor)

  model$parameters[[length(model$parameters)+1]] <- paramInfo
  return(model)
}


# Get AZRmodel parameter information
#
# @param model An AZRmodel
# @param index Index of the parameter in the model
# @return A list with the parameter information
# @examples
# model <- exNovakTyson
# getParameterAZRmodel(model,1)
# @export
getParameterAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("getParameterAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfParametersAZRmodel(model) < index)
    stop("getParameterAZRmodel: value of index larger than number of states.")

  if (index < 1)
    stop("getParameterAZRmodel: value of index should be larger than 0.")

  return(model$parameters[[index]])
}


# Set parameter information in an AZRmodel
#
# Only provided information is updated in the parameter definition.
#
# @param model An AZRmodel
# @param index Numerical index of the parameter to update
# @param name String with parameter name
# @param value Numerical value of parameter
# @param type String with type of SBML element (isSpecie, isParameter, isCompartment)
# @param compartment String with name of SBML compartment for species or outside compartment for compartment)
# @param unittype String with type of unit of SBML species (concentration, amount)
# @param notes String with notes about the parameter
# @param estimate Flag for parameter being estimated
# @param regressor Flag for parameter being a regressor in estimation
# @return An AZRmodel object with updated parameter
# @examples
# model <- exNovakTyson
# setParameterAZRmodel(model,1,value=0.14)
# @export
setParameterAZRmodel <- function(model, index,
                                 name = NULL,
                                 value = NULL,
                                 type = NULL,
                                 compartment = NULL,
                                 unittype = NULL,
                                 notes = NULL,
                                 estimate = FALSE,
                                 regressor = FALSE) {

  if (!is.AZRmodel(model))
    stop("setParameterAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfParametersAZRmodel(model) < index)
    stop("setParameterAZRmodel: value of index larger than number of parameters.")

  if (index < 1)
    stop("setParameterAZRmodel: value of index should be larger than 0.")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("setParameterAZRmodel: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("setParameterAZRmodel: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

  if (!is.null(estimate) && !is.null(regressor) && estimate && regressor)
    stop("addStateAZRmodel: estimate and regressor can not be TRUE at the same time")

  if (!is.null(name))           model$parameters[[index]]$name           <- name
  if (!is.null(value))          model$parameters[[index]]$value          <- value
  if (!is.null(type))           model$parameters[[index]]$type           <- type
  if (!is.null(compartment))    model$parameters[[index]]$compartment    <- compartment
  if (!is.null(unittype))       model$parameters[[index]]$unittype       <- unittype
  if (!is.null(notes))          model$parameters[[index]]$notes          <- notes
  if (!is.null(estimate))       model$parameters[[index]]$estimate       <- estimate
  if (!is.null(regressor))      model$parameters[[index]]$regressor      <- regressor

  return(model)
}


# Delete AZRmodel parameter
#
# @param model An AZRmodel
# @param index Index of the parameter in the model
# @return An AZRmodel with the indexed parameter removed
# @examples
# model <- exNovakTyson
# delParameterAZRmodel(model,1)
# @export
delParameterAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("delParameterAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfParametersAZRmodel(model) < index)
    stop("delParameterAZRmodel: value of index larger than number of parameters.")

  if (index < 1)
    stop("delParameterAZRmodel: value of index should be larger than 0.")

  model$parameters[[index]] <- NULL

  # Need to update the parindex fields in potentially available inputs
  if (getNumberOfInputsAZRmodel(model) > 0) {
    paramnames <- getAllParametersAZRmodel(model)$paramnames
    for (k in 1:getNumberOfInputsAZRmodel(model)) {
      model$inputs[[k]]$parindex = strmatch(paramnames,model$inputs[[k]]$name)
    }
  }

  return(model)
}


###############################################################################
# Variable handling functions
# addVariableAZRmodel
# getVariableAZRmodel
# setVariableAZRmodel
# delVariableAZRmodel
###############################################################################

# Add a new variable to an AZRmodel
#
# @param model An AZRmodel
# @param name String with variable name
# @param formula String with formula for variable
# @param type String with type of SBML element (isSpecie, isVariable, isCompartment)
# @param compartment String with name of SBML compartment for species or outside compartment for compartment)
# @param unittype String with type of unit of SBML species (concentration, amount)
# @param notes String with notes about the variable
# @return An AZRmodel object with appended variable
# @examples
# model <- AZRmodel()
# addVariableAZRmodel(model,'abc','a+b',notes="hello")
# addVariableAZRmodel(model,'F',formula='a+b')
# @export
addVariableAZRmodel <- function(model,
                                 name = NULL,
                                 formula = NULL,
                                 type = NULL,
                                 compartment = NULL,
                                 unittype = NULL,
                                 notes = NULL) {

  if (!is.AZRmodel(model))
    stop("addVariableAZRmodel: model argument is not an AZRmodel")

  if (is.null(name))
    stop("addVariableAZRmodel: name is a required input argument")

  if (is.null(formula))
    stop("addVariableAZRmodel: formula is a required input argument")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("addVariableAZRmodel: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("addVariableAZRmodel: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

  varInfo <- list(name=name, formula=formula,
                   type = type,
                   compartment = compartment,
                   unittype = unittype,
                   notes = notes)

  model$variables[[length(model$variables)+1]] <- varInfo
  return(model)
}


# Get AZRmodel variable information
#
# @param model An AZRmodel
# @param index Index of the variable in the model
# @return A list with the variable information
# @examples
# model <- exNovakTyson
# getVariableAZRmodel(model,1)
# @export
getVariableAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("getVariableAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfVariablesAZRmodel(model) < index)
    stop("getVariableAZRmodel: value of index larger than number of variables.")

  if (index < 1)
    stop("getVariableAZRmodel: value of index should be larger than 0.")

  return(model$variables[[index]])
}


# Set variable information in an AZRmodel
#
# Only provided information is updated in the variable definition.
#
# @param model An AZRmodel
# @param index Numerical index of the variable to update
# @param name String with variable name
# @param formula Numerical value of variable
# @param type String with type of SBML element (isSpecie, isVariable, isCompartment)
# @param compartment String with name of SBML compartment for species or outside compartment for compartment)
# @param unittype String with type of unit of SBML species (concentration, amount)
# @param notes String with notes about the variable
# @return An AZRmodel object with updated variable
# @examples
# model <- exNovakTyson
# model <- setVariableAZRmodel(model,1,formula='Cyclin')
# @export
setVariableAZRmodel <- function(model, index,
                                 name = NULL,
                                 formula = NULL,
                                 type = NULL,
                                 compartment = NULL,
                                 unittype = NULL,
                                 notes = NULL) {

  if (!is.AZRmodel(model))
    stop("setVariableAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfVariablesAZRmodel(model) < index)
    stop("setVariableAZRmodel: value of index larger than number of variables.")

  if (index < 1)
    stop("setVariableAZRmodel: value of index should be larger than 0.")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("setVariableAZRmodel: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("setVariableAZRmodel: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

  if (!is.null(name))           model$variables[[index]]$name           <- name
  if (!is.null(formula))        model$variables[[index]]$formula        <- formula
  if (!is.null(type))           model$variables[[index]]$type           <- type
  if (!is.null(compartment))    model$variables[[index]]$compartment    <- compartment
  if (!is.null(unittype))       model$variables[[index]]$unittype       <- unittype
  if (!is.null(notes))          model$variables[[index]]$notes          <- notes

  return(model)
}


# Delete AZRmodel variable
#
# @param model An AZRmodel
# @param index Index of the variable in the model
# @return An AZRmodel with the indexed variable removed
# @examples
# model <- exNovakTyson
# delVariableAZRmodel(model,1)
# @export
delVariableAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("delVariableAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfVariablesAZRmodel(model) < index)
    stop("delVariableAZRmodel: value of index larger than number of variables.")

  if (index < 1)
    stop("delVariableAZRmodel: value of index should be larger than 0.")

  # Delete variable
  model$variables[[index]] <- NULL

  # Need to update the varindex fields in potentially available outputs
  if (getNumberOfOutputsAZRmodel(model) > 0) {
    varnames <- getAllVariablesAZRmodel(model)$varnames
    for (k in 1:getNumberOfOutputsAZRmodel(model)) {
      model$outputs[[k]]$varindex = strmatch(varnames,model$outputs[[k]]$name)
    }
  }

  return(model)
}


###############################################################################
# Reaction handling functions
# addReactionAZRmodel
# getReactionAZRmodel
# setReactionAZRmodel
# delReactionAZRmodel
###############################################################################

# Add a new reaction to an AZRmodel
#
# @param model An AZRmodel
# @param name String with reaction name
# @param formula String with formula for reaction
# @param notes String with notes about the reaction
# @param reversible Flag for reversible reaction (TRUE) otherwise FALSE
# @param fast Flag for fast reaction (TRUE) otherwise FALSE
# @return An AZRmodel object with appended reaction
# @examples
# model <- AZRmodel()
# addReactionAZRmodel(model,'R1','hello',notes='hello notes',reversible=TRUE)
# addReactionAZRmodel(model,'R2',formula='a+b')
# @export
addReactionAZRmodel <- function(model,
                                name = NULL,
                                formula = NULL,
                                notes = NULL,
                                reversible = FALSE,
                                fast = FALSE) {

  if (!is.AZRmodel(model))
    stop("addReactionAZRmodel: model argument is not an AZRmodel")

  if (is.null(name))
    stop("addReactionAZRmodel: name is a required input argument")

  if (is.null(formula))
    stop("addReactionAZRmodel: formula is a required input argument")

  reacInfo <- list(name=name, formula=formula,
                 notes = notes,
                 reversible = reversible,
                 fast = fast)

  model$reactions[[length(model$reactions)+1]] <- reacInfo
  return(model)
}


# Get AZRmodel reaction information
#
# @param model An AZRmodel
# @param index Index of the reaction in the model
# @return A list with the reaction information
# @examples
# model <- exNovakTyson
# getReactionAZRmodel(model,19)
# @export
getReactionAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("getReactionAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfReactionsAZRmodel(model) < index)
    stop("getReactionAZRmodel: value of index larger than number of reactions.")

  if (index < 1)
    stop("getReactionAZRmodel: value of index should be larger than 0.")

  return(model$reactions[[index]])
}


# Set reaction information in an AZRmodel
#
# Only provided information is updated in the reaction definition.
#
# @param model An AZRmodel
# @param index Numerical index of the reaction to update
# @param name String with reaction name
# @param formula String with formula for reaction
# @param notes String with notes about the reaction
# @param reversible Flag for reversible reaction (TRUE) otherwise FALSE
# @param fast Flag for fast reaction (TRUE) otherwise FALSE
# @return An AZRmodel object with updated reaction
# @examples
# model <- exNovakTyson
# model <- setReactionAZRmodel(model,19,formula='R1')
# @export
setReactionAZRmodel <- function(model, index,
                                name = NULL,
                                formula = NULL,
                                notes = NULL,
                                reversible = NULL,
                                fast = NULL) {

  if (!is.AZRmodel(model))
    stop("setReactionAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfReactionsAZRmodel(model) < index)
    stop("setReactionAZRmodel: value of index larger than number of reactions.")

  if (index < 1)
    stop("setReactionAZRmodel: value of index should be larger than 0.")

  if (!is.null(name))           model$reactions[[index]]$name           <- name
  if (!is.null(formula))        model$reactions[[index]]$formula        <- formula
  if (!is.null(notes))          model$reactions[[index]]$notes          <- notes
  if (!is.null(reversible))     model$reactions[[index]]$reversible     <- reversible
  if (!is.null(fast))           model$reactions[[index]]$fast           <- fast

  return(model)
}


# Delete AZRmodel reaction
#
# @param model An AZRmodel
# @param index Index of the reaction in the model
# @return An AZRmodel with the indexed reaction removed
# @examples
# model <- exNovakTyson
# delReactionAZRmodel(model,1)
# @export
delReactionAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("delReactionAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfReactionsAZRmodel(model) < index)
    stop("delReactionAZRmodel: value of index larger than number of reactions.")

  if (index < 1)
    stop("delReactionAZRmodel: value of index should be larger than 0.")

  model$reactions[[index]] <- NULL

  return(model)
}


###############################################################################
# Function handling functions
# addFunctionAZRmodel
# getFunctionAZRmodel
# setFunctionAZRmodel
# delFunctionAZRmodel
###############################################################################

# Add a new function to an AZRmodel
#
# @param model An AZRmodel
# @param name String with function name
# @param arguments String with comma separated arguments
# @param formula String with formula for function
# @param notes String with notes about the function
# @return An AZRmodel object with appended function
# @examples
# model <- AZRmodel()
# addFunctionAZRmodel(model,'ADD',arguments="x,y",formula="x+y")
# addFunctionAZRmodel(model,'MM',arguments='X,VMAX,KM',formula='VMAX*X/(X+KM)')
# @export
addFunctionAZRmodel <- function(model,
                                name = NULL,
                                arguments = NULL,
                                formula = NULL,
                                notes = NULL) {

  if (!is.AZRmodel(model))
    stop("addFunctionAZRmodel: model argument is not an AZRmodel")

  if (is.null(name))
    stop("addFunctionAZRmodel: name is a required input argument")

  if (is.null(arguments))
    stop("addFunctionAZRmodel: arguments is a required input argument")

  if (is.null(formula))
    stop("addFunctionAZRmodel: formula is a required input argument")

  funInfo <- list(name=name, arguments=arguments, formula=formula, notes = notes)

  model$functions[[length(model$functions)+1]] <- funInfo
  return(model)
}


# Get AZRmodel function information
#
# @param model An AZRmodel
# @param index Index of the function in the model
# @return A list with the function information
# @examples
# model <- AZRmodel()
# model <- addFunctionAZRmodel(model,'ADD',arguments="x,y",formula="x+y")
# getFunctionAZRmodel(model,1)
# @export
getFunctionAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("getFunctionAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfFunctionsAZRmodel(model) < index)
    stop("getFunctionAZRmodel: value of index larger than number of functions.")

  if (index < 1)
    stop("getFunctionAZRmodel: value of index should be larger than 0.")

  return(model$functions[[index]])
}


# Set function information in an AZRmodel
#
# Only provided information is updated in the function definition.
#
# @param model An AZRmodel
# @param index Numerical index of the function to update
# @param name String with function name
# @param arguments String with comma separated arguments
# @param formula String with formula for function
# @param notes String with notes about the function
# @return An AZRmodel object with updated function
# @examples
# model <- AZRmodel()
# model <- addFunctionAZRmodel(model,'ADD',arguments="x,y",formula="x+y")
# setFunctionAZRmodel(model,1,formula='x*y')
# @export
setFunctionAZRmodel <- function(model, index,
                                name = NULL,
                                arguments = NULL,
                                formula = NULL,
                                notes = NULL) {

  if (!is.AZRmodel(model))
    stop("setFunctionAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfFunctionsAZRmodel(model) < index)
    stop("setFunctionAZRmodel: value of index larger than number of functions.")

  if (index < 1)
    stop("setFunctionAZRmodel: value of index should be larger than 0.")

  if (!is.null(name))           model$functions[[index]]$name           <- name
  if (!is.null(arguments))      model$functions[[index]]$arguments      <- arguments
  if (!is.null(formula))        model$functions[[index]]$formula        <- formula
  if (!is.null(notes))          model$functions[[index]]$notes          <- notes

  return(model)
}


# Delete AZRmodel function
#
# @param model An AZRmodel
# @param index Index of the function in the model
# @return An AZRmodel with the indexed function removed
# @examples
# model <- AZRmodel()
# model <- addFunctionAZRmodel(model,'ADD',arguments="x,y",formula="x+y")
# delFunctionAZRmodel(model,1)
# @export
delFunctionAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("delFunctionAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfFunctionsAZRmodel(model) < index)
    stop("delFunctionAZRmodel: value of index larger than number of functions.")

  if (index < 1)
    stop("delFunctionAZRmodel: value of index should be larger than 0.")

  model$functions[[index]] <- NULL

  return(model)
}


###############################################################################
# Algebraic handling functions
# addAlgebraicAZRmodel
# getAlgebraicAZRmodel
# setAlgebraicAZRmodel
# delAlgebraicAZRmodel
###############################################################################

# Add a new algebraic state to an AZRmodel
#
# @param model An AZRmodel
# @param name String with algebraic state name
# @param IC Initial condition (numeric or string)
# @param formula String with formula for algebraic state (evaluated to 0)
# @param type String with type of SBML element (isSpecie, isParameter, isCompartment)
# @param compartment String with name of SBML compartment for species or outside compartment for compartment)
# @param unittype String with type of unit (concentration, amount)
# @param notes String with notes about the algebraic state
# @return An AZRmodel object with appended algebraic state variable
# @examples
# model <- AZRmodel()
# model <- addAlgebraicAZRmodel(model,'X',IC=2,'A+B+X')
# model <- addAlgebraicAZRmodel(model,formula="A+B+X")
# @export
addAlgebraicAZRmodel <- function(model,
                             name = NULL,
                             IC = NULL,
                             formula = NULL,
                             type = NULL,
                             compartment = NULL,
                             unittype = NULL,
                             notes = NULL) {

  if (!is.AZRmodel(model))
    stop("addAlgebraicAZRmodel: model argument is not an AZRmodel")

  if (is.null(name) && !is.null(IC))
    stop("addAlgebraicAZRmodel: name is a required input argument if an IC is defined")

  if (is.null(IC) && !is.null(name))
    stop("addAlgebraicAZRmodel: IC is a required input argument if a name is defined")

  if (is.null(formula))
    stop("addAlgebraicAZRmodel: formula is a required input argument")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("addAlgebraicAZRmodel: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("addAlgebraicAZRmodel: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

  algebraicInfo <- list(name=name, IC=IC, formula=formula,
                   type = type,
                   compartment = compartment,
                   unittype = unittype,
                   notes = notes)

  model$algebraic[[length(model$algebraic)+1]] <- algebraicInfo
  return(model)
}


# Get AZRmodel algebraic state information
#
# @param model An AZRmodel
# @param index Index of the algebraic state in the model
# @return A list with the algebraic state information
# @examples
# model <- AZRmodel()
# model <- addAlgebraicAZRmodel(model,'X',IC=2,'A+B+X')
# getAlgebraicAZRmodel(model,1)
# @export
getAlgebraicAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("getAlgebraicAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfAlgebraicAZRmodel(model) < index)
    stop("getAlgebraicAZRmodel: value of index larger than number of algebraic states.")

  if (index < 1)
    stop("getAlgebraicAZRmodel: value of index should be larger than 0.")

  return(model$algebraic[[index]])
}


# Set algebraic state information in an AZRmodel
#
# Only provided information is updated in the algebraic state definition.
#
# @param model An AZRmodel
# @param index Numerical index of the algebraic state to update
# @param name String with algebraic state name
# @param IC Initial condition (numeric or string)
# @param formula String with formula for algebraic state (evaluated to 0)
# @param type String with type of SBML element (isSpecie, isParameter, isCompartment)
# @param compartment String with name of SBML compartment for species or outside compartment for compartment)
# @param unittype String with type of unit (concentration, amount)
# @param notes String with notes about the algebraic state
# @return An AZRmodel object with updated algebraic state variable
# @examples
# model <- AZRmodel()
# model <- addAlgebraicAZRmodel(model,'X',IC=2,'A+B+X')
# setAlgebraicAZRmodel(model,1,IC=0.14)
# @export
setAlgebraicAZRmodel <- function(model, index,
                             name = NULL,
                             IC = NULL,
                             formula = NULL,
                             type = NULL,
                             compartment = NULL,
                             unittype = NULL,
                             notes = NULL) {

  if (!is.AZRmodel(model))
    stop("setAlgebraicAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfAlgebraicAZRmodel(model) < index)
    stop("setAlgebraicAZRmodel: value of index larger than number of algebraic states.")

  if (index < 1)
    stop("setAlgebraicAZRmodel: value of index should be larger than 0.")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("setAlgebraicAZRmodel: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("setAlgebraicAZRmodel: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

  if (!is.null(name))           model$algebraic[[index]]$name           <- name
  if (!is.null(IC))             model$algebraic[[index]]$IC             <- IC
  if (!is.null(formula))        model$algebraic[[index]]$formula        <- formula
  if (!is.null(type))           model$algebraic[[index]]$type           <- type
  if (!is.null(compartment))    model$algebraic[[index]]$compartment    <- compartment
  if (!is.null(unittype))       model$algebraic[[index]]$unittype       <- unittype
  if (!is.null(notes))          model$algebraic[[index]]$notes          <- notes

  return(model)
}


# Delete AZRmodel algebraic state
#
# @param model An AZRmodel
# @param index Index of the algebraic state in the model
# @return An AZRmodel with the indexed algebraic state removed
# @examples
# model <- AZRmodel()
# model <- addAlgebraicAZRmodel(model,'X',IC=2,'A+B+X')
# delAlgebraicAZRmodel(model,1)
# @export
delAlgebraicAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("delAlgebraicAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfAlgebraicAZRmodel(model) < index)
    stop("delAlgebraicAZRmodel: value of index larger than number of algebraic states.")

  if (index < 1)
    stop("delAlgebraicAZRmodel: value of index should be larger than 0.")

  model$algebraic[[index]] <- NULL

  return(model)
}


###############################################################################
# Event handling functions
# addEventAZRmodel
# getEventAZRmodel
# setEventAZRmodel
# delEventAZRmodel
###############################################################################

# Add a new event to an AZRmodel
#
# This function only adds an event but not event assignments. For these separate
# handling functions are available.
#
# @param model An AZRmodel
# @param name String with event name
# @param trigger String with event trigger information
# @param notes String with notes about the event
# @return An AZRmodel object with appended event
# @examples
# model <- AZRmodel()
# addEventAZRmodel(model,'event1',trigger="lt(X,2)")
# addEventAZRmodel(model,'event2',trigger="gt(time,10)")
# @export
addEventAZRmodel <- function(model,
                                 name = NULL,
                                 trigger = NULL,
                                 notes = NULL) {

  if (!is.AZRmodel(model))
    stop("addEventAZRmodel: model argument is not an AZRmodel")

  if (is.null(name))
    stop("addEventAZRmodel: name is a required input argument")

  if (is.null(trigger))
    stop("addEventAZRmodel: trigger is a required input argument")

  eventInfo <- list(name=name, trigger=trigger, assignment=NULL, notes=notes)

  model$events[[length(model$events)+1]] <- eventInfo
  return(model)
}


# Get AZRmodel event information
#
# @param model An AZRmodel
# @param index Index of the event in the model
# @return A list with the event information
# @examples
# model <- AZRmodel()
# model <- addEventAZRmodel(model,'event1',trigger="lt(X,2)")
# getEventAZRmodel(model,1)
# @export
getEventAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("getEventAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfEventsAZRmodel(model) < index)
    stop("getEventAZRmodel: value of index larger than number of events.")

  if (index < 1)
    stop("getEventAZRmodel: value of index should be larger than 0.")

  return(model$events[[index]])
}


# Set event information in an AZRmodel
#
# Only provided information is updated in the event definition.
#
# @param model An AZRmodel
# @param index Numerical index of the event to update
# @param name String with event name
# @param trigger String with event trigger information
# @param notes String with notes about the event
# @return An AZRmodel object with updated event variable
# @examples
# model <- AZRmodel()
# model <- addEventAZRmodel(model,'event1',trigger="lt(X,2)")
# setEventAZRmodel(model,1,name='myEvent')
# @export
setEventAZRmodel <- function(model, index,
                                 name = NULL,
                                 trigger = NULL,
                                 notes = NULL) {

  if (!is.AZRmodel(model))
    stop("setEventAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfEventsAZRmodel(model) < index)
    stop("setEventAZRmodel: value of index larger than number of events.")

  if (index < 1)
    stop("setEventAZRmodel: value of index should be larger than 0.")

  if (!is.null(name))           model$events[[index]]$name           <- name
  if (!is.null(trigger))        model$events[[index]]$trigger        <- trigger
  if (!is.null(notes))          model$events[[index]]$notes          <- notes

  return(model)
}


# Delete AZRmodel event
#
# @param model An AZRmodel
# @param index Index of the event in the model
# @return An AZRmodel with the indexed event removed
# @examples
# model <- AZRmodel()
# model <- addEventAZRmodel(model,'event1',trigger="lt(X,2)")
# delEventAZRmodel(model,1)
# @export
delEventAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("delEventAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfEventsAZRmodel(model) < index)
    stop("delEventAZRmodel: value of index larger than number of events.")

  if (index < 1)
    stop("delEventAZRmodel: value of index should be larger than 0.")

  model$events[[index]] <- NULL

  return(model)
}


###############################################################################
# Event Assignment handling functions
# addEventAssignmentAZRmodel
# getEventAssignmentAZRmodel
# setEventAssignmentAZRmodel
# delEventAssignmentAZRmodel
###############################################################################

# Add a new event assignment to an event in an AZRmodel
#
# This function adds an event assignment to a specified event. It does not create
# an event.
#
# @param model An AZRmodel
# @param eventindex Event index in the AZRmodel
# @param variable String with name of variable to change by the event
# @param formula String with event assignment formula
# @return An AZRmodel object with appended event assignment to the specified event
# @examples
# model <- AZRmodel()
# model <- addEventAZRmodel(model,'event1',trigger="lt(X,2)")
# model <- addEventAssignmentAZRmodel(model,eventindex=1,'A','A+5')
# @export
addEventAssignmentAZRmodel <- function(model,
                                       eventindex = NULL,
                                       variable = NULL,
                                       formula = NULL) {

  if (!is.AZRmodel(model))
    stop("addEventAssignmentAZRmodel: model argument is not an AZRmodel")

  if (is.null(eventindex))
    stop("addEventAssignmentAZRmodel: eventindex is a required input argument")

  if (is.null(variable))
    stop("addEventAssignmentAZRmodel: variable is a required input argument")

  if (is.null(formula))
    stop("addEventAssignmentAZRmodel: formula is a required input argument")

  eventAssignmentInfo <- list(variable=variable, formula=formula)

  model$events[[eventindex]]$assignment[[length(model$events[[eventindex]]$assignment)+1]] <- eventAssignmentInfo
  return(model)
}


# Get AZRmodel event assignment information
#
# @param model An AZRmodel
# @param eventindex Index of the event in the model
# @param index Index of the event assignment to return
# @return A list with the event information
# @examples
# model <- AZRmodel()
# model <- addEventAZRmodel(model,'event1',trigger="lt(X,2)")
# model <- addEventAssignmentAZRmodel(model,eventindex=1,'A','A+5')
# getEventAssignmentAZRmodel(model,1,1)
# @export
getEventAssignmentAZRmodel <- function(model, eventindex, index) {

  if (!is.AZRmodel(model))
    stop("getEventAssignmentAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfEventsAZRmodel(model) < eventindex)
    stop("getEventAssignmentAZRmodel: value of eventindex larger than number of events.")

  if (eventindex < 1)
    stop("getEventAssignmentAZRmodel: value of eventindex should be larger than 0.")

  if (getNumberOfEventassignmentsAZRmodel(model,eventindex) < index)
    stop("getEventAssignmentAZRmodel: value of index larger than number of event assignments in selected event.")

  if (index < 1)
    stop("getEventAssignmentAZRmodel: value of index should be larger than 0.")

  return(model$events[[eventindex]]$assignment[[index]])
}


# Set event assignment information in an AZRmodel
#
# Only provided information is updated in the event assignment definition.
#
# @param model An AZRmodel
# @param eventindex Numerical index of the event to update
# @param index Numerical index of the event assignment to update
# @param variable String with name of variable to change by the event
# @param formula String with event assignment formula
# @return An AZRmodel object with updated event assignment
# @examples
# model <- AZRmodel()
# model <- addEventAZRmodel(model,'event1',trigger="lt(X,2)")
# model <- addEventAssignmentAZRmodel(model,eventindex=1,'A','A+5')
# setEventAssignmentAZRmodel(model,1,1,variable='b')
# @export
setEventAssignmentAZRmodel <- function(model, eventindex, index,
                             variable = NULL,
                             formula = NULL) {

  if (!is.AZRmodel(model))
    stop("setEventAssignmentAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfEventsAZRmodel(model) < eventindex)
    stop("setEventAssignmentAZRmodel: value of eventindex larger than number of events.")

  if (eventindex < 1)
    stop("setEventAssignmentAZRmodel: value of eventindex should be larger than 0.")

  if (getNumberOfEventassignmentsAZRmodel(model,eventindex) < index)
    stop("setEventAssignmentAZRmodel: value of index larger than number of event assignments in selected event.")

  if (index < 1)
    stop("setEventAssignmentAZRmodel: value of index should be larger than 0.")

  if (!is.null(variable))       model$events[[eventindex]]$assignment[[index]]$variable       <- variable
  if (!is.null(formula))        model$events[[eventindex]]$assignment[[index]]$formula        <- formula

  return(model)
}


# Delete AZRmodel event assignment
#
# @param model An AZRmodel
# @param eventindex Numerical index of the event to delete
# @param index Numerical index of the event assignment to delete
# @return An AZRmodel with the indexed event assignment removed
# @examples
# model <- AZRmodel()
# model <- addEventAZRmodel(model,'event1',trigger="lt(X,2)")
# model <- addEventAssignmentAZRmodel(model,eventindex=1,'A','A+5')
# delEventAssignmentAZRmodel(model,1,1)
# @export
delEventAssignmentAZRmodel <- function(model, eventindex, index) {

  if (!is.AZRmodel(model))
    stop("delEventAssignmentAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfEventsAZRmodel(model) < eventindex)
    stop("delEventAssignmentAZRmodel: value of eventindex larger than number of events.")

  if (eventindex < 1)
    stop("delEventAssignmentAZRmodel: value of eventindex should be larger than 0.")

  if (getNumberOfEventassignmentsAZRmodel(model,eventindex) < index)
    stop("delEventAssignmentAZRmodel: value of index larger than number of event assignments in selected event.")

  if (index < 1)
    stop("delEventAssignmentAZRmodel: value of index should be larger than 0.")

  model$events[[eventindex]]$assignment[[index]] <- NULL

  return(model)
}


###############################################################################
# Input handling functions
# addInputAZRmodel
# getInputAZRmodel
# delInputAZRmodel
###############################################################################

# Add a new input to an AZRmodel
#
# This function adds an input to an AZRmodel. No names can be specified.
# First input added is INPUT1, second INPUT2, etc. Inputs can only be added to
# ODEs and will have an impact on the rate of change of the specified
# state variable. An input can be added on several states. For each state to
# which an input is added a factor needs to be defined in the format "+factor".
# A "factor" can consist of a mathematical expression of numerical values and
# parameters, defined in the AZRmodel. Adding an input will add a parameter
# to the AZRmodel with the same name as the input and initialize it to 0.
# It will also add the required input terms into the ODEs.
#
# @param model An AZRmodel
# @param factors  String or vector of strings with factors for the input. It is important
#                 that factors do not consist of additive terms. Allowed: factors=c("(1-F)","F").
#                 Not allowed: factors=c("1-F",F). It is allowed to have a "+" as first character
#                 in a factor. Example: factors=c("+(1-F)","+F"). If the "+" is not present, it will
#                 be added when storing the information in the model.
# @param stateindex Scalar or vector with index/indices of states on which to add the input
# @return An AZRmodel object with appended input
# @examples
# model <- AZRmodel()
# model <- addStateAZRmodel(model,'Cyclin',IC=0.3,'Re1-Re2')
# model <- addStateAZRmodel(model,'Cyclin2',IC=0.3,'Re2-Re3')
# model <- addInputAZRmodel(model,stateindex=2)
# model <- addInputAZRmodel(model,factors='F',stateindex=2)
# model <- addInputAZRmodel(model,factors=c('+F','(1-F)'),stateindex=c(1,2))
# @export
addInputAZRmodel <- function(model,
                             stateindex = NULL,
                             factors = "+1") {

  if (!is.AZRmodel(model))
    stop("addInputAZRmodel: model argument is not an AZRmodel")

  if (is.null(stateindex))
    stop("addInputAZRmodel: stateindex is a required input argument")

  if (length(factors) != length(stateindex))
    stop("addInputAZRmodel: stateindex and factors need to have same number of elements")

  if (max(stateindex) > getNumberOfStatesAZRmodel(model))
    stop("addInputAZRmodel: values of stateindices exceed number of states in the model")

  if (min(stateindex) < 1)
    stop("addInputAZRmodel: values of stateindices need to be larger than 0")

  name <- paste("INPUT", getNumberOfInputsAZRmodel(model)+1, sep="")
  for (k in 1:length(factors)) {
    factors[k] <- strtrim(factors[k])
    if (substr(factors[k],1,1)!="+") factors[k] <- paste("+",factors[k],sep="")
  }

  terms <- factors
  for (k in 1:length(terms)) terms[k] <- paste(strremWhite(factors[k]),"*",name, sep="")
  model <- addParameterAZRmodel(model,name=name,value=0)
  parindex <- getNumberOfParametersAZRmodel(model)
  for (k in 1:length(stateindex))
    model$states[[stateindex[k]]]$ODE <- paste(model$states[[stateindex[k]]]$ODE,terms[k],sep="")

  inputInfo <- list(name=name, factors=factors, terms=terms, stateindex=stateindex, parindex=parindex)

  model$inputs[[getNumberOfInputsAZRmodel(model)+1]] <- inputInfo
  return(model)
}


# Get AZRmodel input information
#
# @param model An AZRmodel
# @param index Index of the input in the model
# @return A list with the input information
# @examples
# model <- AZRmodel()
# model <- addStateAZRmodel(model,'Cyclin',IC=0.3,'Re1-Re2')
# model <- addStateAZRmodel(model,'Cyclin2',IC=0.3,'Re2-Re3')
# model <- addInputAZRmodel(model,stateindex=2)
# getInputAZRmodel(model,1)
# @export
getInputAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("getInputAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfInputsAZRmodel(model) < index)
    stop("getInputAZRmodel: value of index larger than number of inputs.")

  if (index < 1)
    stop("getInputAZRmodel: value of index should be larger than 0.")

  return(model$input[[index]])
}

# Delete AZRmodel input
#
# @param model An AZRmodel
# @param index Numerical index of the input to delete
# @return An AZRmodel with the indexed input removed (parameter will be removed as well)
# @examples
# model <- AZRmodel()
# model <- addStateAZRmodel(model,'Cyclin',IC=0.3,'Re1-Re2')
# model <- addStateAZRmodel(model,'Cyclin2',IC=0.3,'Re2-Re3')
# model <- addInputAZRmodel(model,stateindex=2)
# delInputAZRmodel(model,1)
# @export
delInputAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("delInputAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfInputsAZRmodel(model) < index)
    stop("delInputAZRmodel: value of index larger than number of inputs.")

  if (index < 1)
    stop("delInputAZRmodel: value of index should be larger than 0.")

  # Remove INPUT term in ODE
  si <- model$inputs[[index]]$stateindex
  for (k in 1:length(si))
    model$states[[si[k]]]$ODE <- strrep(strremWhite(model$states[[si[k]]]$ODE),model$inputs[[index]]$terms[[k]],"")

  # Save parindex of input
  parindex <- model$inputs[[index]]$parindex

  # Remove input
  model$inputs[[index]] <- NULL

  # Remove INPUT parameter
  model <- delParameterAZRmodel(model,parindex)

  return(model)
}

###############################################################################
# Output handling functions
# addOutputAZRmodel
# getOutputAZRmodel
# delOutputAZRmodel
###############################################################################

# Add a new output to an AZRmodel
#
# This function adds an output to an AZRmodel. No names can be specified.
# First output added is OUTPUT1, second OUTPUT2, etc. Outputs are added as
# variables and their formula can contain expressions containing states, parameters,
# and already defined variables.
#
# @param model An AZRmodel
# @param formula  Formula for the created output variable
# @param notes String with notes for the output variable
# @return An AZRmodel object with appended input
# @examples
# model <- AZRmodel()
# addOutputAZRmodel(model,formula='Ac/Vc')
# @export
addOutputAZRmodel <- function(model,
                              formula = NULL,
                              notes = NULL) {

  if (!is.AZRmodel(model))
    stop("addOutputAZRmodel: model argument is not an AZRmodel")

  if (is.null(formula))
    stop("addOutputAZRmodel: formula is a required input argument")

  name <- paste("OUTPUT", getNumberOfOutputsAZRmodel(model)+1, sep="")

  model <- addVariableAZRmodel(model,name,formula,notes=notes)

  varindex <- getNumberOfVariablesAZRmodel(model)

  outputInfo <- list(name=name, formula=formula, notes=notes, varindex=varindex)

  model$outputs[[getNumberOfOutputsAZRmodel(model)+1]] <- outputInfo
  return(model)
}


# Get AZRmodel output information
#
# @param model An AZRmodel
# @param index Index of the output in the model
# @return A list with the output information
# @examples
# model <- AZRmodel()
# model <- addOutputAZRmodel(model,formula='Ac/Vc')
# getOutputAZRmodel(model,1)
# @export
getOutputAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("getOutputAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfOutputsAZRmodel(model) < index)
    stop("getOutputAZRmodel: value of index larger than number of outputs.")

  if (index < 1)
    stop("getOutputAZRmodel: value of index should be larger than 0.")

  return(model$output[[index]])
}

# Delete AZRmodel output
#
# @param model An AZRmodel
# @param index Numerical index of the output to delete
# @return An AZRmodel with the indexed output removed (variable will be removed as well)
# @examples
# model <- AZRmodel()
# model <- addOutputAZRmodel(model,formula='Ac/Vc')
# delOutputAZRmodel(model,1)
# @export
delOutputAZRmodel <- function(model, index) {

  if (!is.AZRmodel(model))
    stop("delOutputAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfOutputsAZRmodel(model) < index)
    stop("delOutputAZRmodel: value of index larger than number of outputs.")

  if (index < 1)
    stop("delOutputAZRmodel: value of index should be larger than 0.")

  # Get variable index related with output
  varindex = model$outputs[[index]]$varindex

  # First delete output
  model$outputs[[index]] <- NULL

  # Then delete variable related to the output
  model <- delVariableAZRmodel(model,varindex)

  return(model)
}


###############################################################################
# Get all element functions
# getAllStatesAZRmodel
# getAllParametersAZRmodel
# getAllVariablesAZRmodel
# getAllReactionsAZRmodel
###############################################################################

# Get information about all states in an AZRmodel
#
# @param model An AZRmodel
# @return A list with the following components:
# \item{statenames}{Vector with state names}
# \item{stateICs}{Vector with initial conditions}
# \item{stateODEs}{Vector with ODEs}
# @examples
# model <- exNovakTyson
# getAllStatesAZRmodel(model)
# @export
getAllStatesAZRmodel <- function(model) {

  if (!is.AZRmodel(model))
    stop("getAllStatesAZRmodel: model argument is not an AZRmodel")

  statenames <- c()
  stateICs <- c()
  stateODEs <- c()
  if (getNumberOfStatesAZRmodel(model) > 0) {
    for (k in 1:getNumberOfStatesAZRmodel(model)) {
      x <- getStateAZRmodel(model,k)
      statenames[k] <- x$name
      stateICs[k] <- x$IC
      stateODEs[k] <- x$ODE
    }
  }
  names(statenames) <- statenames
  names(stateICs) <- statenames
  names(stateODEs) <- statenames

  return(list(statenames=statenames,stateICs=stateICs,stateODEs=stateODEs))
}

# Get information about all parameters in an AZRmodel
#
# @param model An AZRmodel
# @return A list with the following components:
# \item{paramnames}{Vector with parameter names}
# \item{paramvalues}{Vector with parameter values}
# \item{paramestimate}{Vector with parameter estimation flags}
# \item{paramregressor}{Vector with parameter regressor flags}
# @examples
# model <- exNovakTyson
# getAllParametersAZRmodel(model)
# @export
getAllParametersAZRmodel <- function(model) {

  if (!is.AZRmodel(model))
    stop("getAllParametersAZRmodel: model argument is not an AZRmodel")

  paramnames <- c()
  paramvalues <- c()
  paramestimate <- c()
  paramregressor <- c()
  if (getNumberOfParametersAZRmodel(model) > 0) {
    for (k in 1:getNumberOfParametersAZRmodel(model)) {
      x <- getParameterAZRmodel(model,k)
      paramnames[k] <- x$name
      paramvalues[k] <- x$value
      paramestimate[k] <- x$estimate
      paramregressor[k] <- x$regressor
    }
  }
  names(paramnames) <- paramnames
  names(paramvalues) <- paramnames
  names(paramestimate) <- paramnames
  names(paramregressor) <- paramnames

  return(list(paramnames=paramnames,paramvalues=paramvalues,paramestimate=paramestimate,paramregressor=paramregressor))
}

# Get information about all variables in an AZRmodel
#
# @param model An AZRmodel
# @return A list with the following components:
# \item{varnames}{Vector with variable names}
# \item{varformulas}{Vector with variable formulas}
# @examples
# model <- exNovakTyson
# getAllVariablesAZRmodel(model)
# @export
getAllVariablesAZRmodel <- function(model) {

  if (!is.AZRmodel(model))
    stop("getAllVariablesAZRmodel: model argument is not an AZRmodel")

  varnames <- c()
  varformulas <- c()
  if (getNumberOfVariablesAZRmodel(model) > 0) {
    for (k in 1:getNumberOfVariablesAZRmodel(model)) {
      x <- getVariableAZRmodel(model,k)
      varnames[k] <- x$name
      varformulas[k] <- x$formula
    }
  }
  names(varnames) <- varnames
  names(varformulas) <- varnames

  return(list(varnames=varnames,varformulas=varformulas))
}

# Get information about all reactions in an AZRmodel
#
# @param model An AZRmodel
# @return A list with the following components:
# \item{reacnames}{Vector with reaction names}
# \item{reacformulas}{Vector with reaction formulas}
# \item{reacfast}{Flag for fast reactions}
# \item{reacreversible}{Flag for reversible reactions}
# @examples
# model <- exNovakTyson
# getAllReactionsAZRmodel(model)
# @export
getAllReactionsAZRmodel <- function(model) {

  if (!is.AZRmodel(model))
    stop("getAllReactionsAZRmodel: model argument is not an AZRmodel")

  reacnames <- c()
  reacformulas <- c()
  reacfast <- c()
  reacreversible <- c()
  if (getNumberOfReactionsAZRmodel(model) > 0) {
    for (k in 1:getNumberOfReactionsAZRmodel(model)) {
      x <- getReactionAZRmodel(model,k)
      reacnames[k] <- x$name
      reacformulas[k] <- x$formula
      reacfast[k] <- x$fast
      reacreversible[k] <- x$reversible
    }
  }
  names(reacnames) <- reacnames
  names(reacformulas) <- reacnames
  names(reacfast) <- reacnames
  names(reacreversible) <- reacnames
  return(list(reacnames=reacnames,reacformulas=reacformulas,reacfast=reacfast,reacreversible=reacreversible))
}


###############################################################################
# Rename model elements
###############################################################################
# Allows to rename elements in the whole model
#
# This function is only allowed to be run when NO simulation functions are attached
# yet. This is checked and an error is returned otherwise.
#
# @param model AZRmodel
# @param origStrings vector of original Strings
# @param newStrings vector or new strings
# @return model updated model

renameElementsAZRmodel <- function(model, origStrings, newStrings) {

  if (!is.null(attr(model,"ODEsim")))
    stop("renameElementsAZRmodel: simulation functions attached - not allowed.")

  # Get temporary text file name
  tempfilename <- paste(tempfile(),".txt",sep="")
  # Export model to temporary text file - use low level functions to avoid issues
  # with simfunction handling !!!
  exportTxtAZRmodel(model,filename=tempfilename)
  # Load text file
  content <- fileread(tempfilename)
  # Exchange strings
  searchStrings <- paste("\\b",origStrings,"\\b",sep="")
  for (k in 1:length(searchStrings)) {
    content <- gsub(searchStrings[k], newStrings[k], content)
  }
  # Save modified model
  filewrite(content,tempfilename)
  # Load model (without generation of simulation functions)
  model <- importTxtAZRmodel(AZRmodel(),tempfilename)
  # Delete temp file
  unlink(tempfilename)
  # Return model
  return(model)
}


###############################################################################
# Replace text in AZRmodel
###############################################################################
# Replace text in AZRmodel
#
# This function can be used to replace a text piece in an AZRmodel.
# Useful to remove "AZRsim:::" strings before export of a model.
# No check will be done for a whole word - so be careful!
#
# @param model AZRmodel
# @param origString string to replace
# @param newString new string
# @return model updated model

replaceTextAZRmodel <- function(model, origString, newString) {

  # Check if states present in the model
  if (getNumberOfStatesAZRmodel(model) == 0)
    return(model)

  # Get temporary text file name
  tempfilename <- paste(tempfile(),".txt",sep="")
  # Export model to temporary text file
  exportTxtAZRmodel(model,filename=tempfilename)
  # Load text file
  content <- fileread(tempfilename)
  # Exchange string
  content <- gsub(origString, newString, content, fixed=TRUE)
  # Save modified model
  filewrite(content,tempfilename)
  # Load model (without generation of simulation functions)
  model <- importTxtAZRmodel(AZRmodel(),tempfilename)
  # Delete temp file
  unlink(tempfilename)
  # Return model
  return(model)
}
