###############################################################################
###############################################################################
# This file contains the create_model function and simple methods for handling azrmod objects
###############################################################################
###############################################################################


###############################################################################
# azrmod: function returing S3 object "azrmod"
###############################################################################
#' azrmod class
#'
#' This function creates an object of class \code{azrmod} which represents dynamic,
#' ODE based models that can be simulated. The class constructor function
#' \code{azrmod} creates such objects from specifc text file descriptions of
#' these models. If no input argument is provided, an empty azrmod object is returned.
#'
#' @param input A string with the full filename (including absolute or
#'   relative path) to a text file, describing an azrmod.
#' @param simFlag If set to \code{TRUE} (default) then the model will be prepared for
#'   simulations (dosing, constraints, etc.) and all required sim functions will
#'   be generated and added to the azrmod attributes. Choose this option when
#'   you want to simulate. In cases were the model structure and the information
#'   about inputs is important, then choose \code{FALSE}.
#' @return An object of class azr_model. If input is NULL, an empty object is returned.
#' @examples
#' empty_mod <- create_model()
#' filename <- system.file(package="AZRsim","examples","NovakTyson.txt")
#' mod <- create_model(filename)
#' @export

create_model <- function (input=NULL,simFlag=TRUE) {

  #################################
  # Handle optional input arguments
  #################################
  # No optional input arguments present

  #################################
  # Checking of input arguments
  #################################
  if(!is.null(input) && !file.exists(input)) {
    # input does not point to an existing file on the file system => error
    stop("Provided path argument does not point to a file on the filesystem.")
  }

  #################################
  # Initialize empty AZRmodel
  #################################
  model <- azrmod_template()

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
  check_azrmod(model)

  #################################
  # Add the original model as attribute (before generation of simulation functions
  # and required update of the model
  #################################
  attr(model,"originalModel") <- model

  #################################
  # Update the models math to for simulation purposes
  # Constraints
  # Inputs / Dosing
  # Generation of C code model and compilation
  #################################
  if (simFlag) {
    model <- genSimFunctions(model)
  }

  #################################
  # Return model
  #################################
  # construct the model object
  class(model) <- "azrmod"
  return(model)
}


###############################################################################
# export_azrmod: exports an AZRmodel as .txt or .txtbc file
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
#' export_azrmod(model)
#' filename <- system.file(package="AZRsim","examples","NovakTyson.txt")
#' model <- AZRmodel(filename)
#' export_azrmod(model,"filename")
#' export_azrmod(model,'filename',useBC=TRUE)
#' @export

export_azrmod <- function (model, filename=NULL, useBC=FALSE, useSIM=FALSE) {

  if (!is_azrmod(model))
    stop("export_azrmod: input argument is not an AZRmodel")

  if (!useSIM && !is.null(attr(model,"originalModel"))) model <- attr(model,"originalModel")

  # Convert to either TEXT or TEXTBC
  if (!useBC)
    exportTxtAZRmodel(model,filename)
  else
    exportTxtBcAZRmodel(model,filename)
}


###############################################################################
# azrmod_template
###############################################################################
# Creates an empty AZRmodel object
#
# @return An empty AZRmodel object.
# @examples
# azrmod_template()

azrmod_template <- function () {

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
  class(model) <- "azrmod"
  return(model)
}


###############################################################################
# Generic function overload (is) - need to be exported
###############################################################################
#' Check function if it is an azrmod
#'
#' @param input azrmod object
#' @return TRUE or FALSE
#' @examples
#'   is_azrmod(AZRmodel())
#' @export
is_azrmod <- function(input) {
  methods::is(input,"azrmod")
}

###############################################################################
# Generic function overload (print) - need to be exported
###############################################################################
#' Generic print function for AZRmodels
#'
#' @param x azrmod object
#' @param ... Additional unused arguments
#' @examples
#'   print(azrmod())
#' @export
print.azrmod <- function(x, ...) {
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
  if (!has_only_numeric_ic(x))
    cat("\tNon-numeric initial conditions are present in the model.\n")
  if (has_fast_reactions(x))
    cat("\tFast reactions defined in the model. Such reactions can not be
        handled by other functions. Please check if these reactions
        are really required. If you need them, contact us (info@intiquan.com
        and we will see what we can do! In the meantime, consider the use
        of IQM Tools (http://www.intiquan.com/iqm-tools/)\n")
}


###############################################################################
# check_azrmod
###############################################################################
# Checks names of model components and issues warnings or errors
check_azrmod <- function(model) {
  stateNames <- toupper(getAllStatesAZRmodel(model)$statenames)
  paramNames <- toupper(getAllParametersAZRmodel(model)$paramnames)
  varNames   <- toupper(getAllVariablesAZRmodel(model)$varnames)
  reacNames  <- toupper(getAllReactionsAZRmodel(model)$reacnames)

  # Check for single char name (should be avoided)
  if (sum(as.numeric(nchar(stateNames) == 1)) > 0)
    warning("check_azrmod: AZRmodel contains state names with a single character name. Try to avoid that if you plan to use NONMEM or MONOLIX")
  if (sum(as.numeric(nchar(paramNames) == 1)) > 0)
    warning("check_azrmod: AZRmodel contains parameter names with a single character name. Try to avoid that if you plan to use NONMEM or MONOLIX")
  if (sum(as.numeric(nchar(varNames) == 1)) > 0)
    warning("check_azrmod: AZRmodel contains variable names with a single character name. Try to avoid that if you plan to use NONMEM or MONOLIX")
  if (sum(as.numeric(nchar(reacNames) == 1)) > 0)
    warning("check_azrmod: AZRmodel contains reaction names with a single character name. Try to avoid that if you plan to use NONMEM or MONOLIX")

  # Define reserved words - case insensitive
  reservedWords <- toupper(c("F","G","H","gt","ge","lt","le","mod","and","or","multiply",
                             "piecewise","interp0","interp1","interpcs", "default", "F1", "F2"))

  # Check if model contains elements with names matching reserved words
  if (length(intersect(stateNames,reservedWords)) > 0)
    stop(paste("check_azrmod: model contains the following state name(s) that is(are) reserved word(s): ",intersect(stateNames,reservedWords),sep=""))
  if (length(intersect(paramNames,reservedWords)) > 0)
    stop(paste("check_azrmod: model contains the following parameter name(s) that is(are) reserved word(s): ",intersect(paramNames,reservedWords),sep=""))
  if (length(intersect(varNames,reservedWords)) > 0)
    stop(paste("check_azrmod: model contains the following variable name(s) that is(are) reserved word(s): ",intersect(varNames,reservedWords),sep=""))
  if (length(intersect(reacNames,reservedWords)) > 0)
    stop(paste("check_azrmod: model contains the following reaction name(s) that is(are) reserved word(s): ",intersect(reacNames,reservedWords),sep=""))
}


###############################################################################
# has_only_numeric_ic
###############################################################################
# Checks if the model contains only numeric initial conditions
#
# @param model An AZRmodel object
# @return TRUE if model contains non-numerical initial conditions, FALSE otherwise
# @examples
# model <- exNovakTyson
# has_only_numeric_ic(model)
has_only_numeric_ic <- function(model) {

  if (!is_azrmod(model))
    stop("has_only_numeric_ic: input argument is not an AZRmodel")

  if (getNumberOfStatesAZRmodel(model) < 1) return(TRUE)

  for (k in 1:getNumberOfStatesAZRmodel(model))
    if (!isnumericVector(model$state[[k]]$IC))
      return(FALSE)

  return(TRUE)
}

###############################################################################
# has_fast_reactions
###############################################################################
# Checks if the model contains fast reactions
#
# @param model An AZRmodel object
# @return TRUE if model contains fast reactions, FALSE otherwise
# @examples
# model <- exNovakTyson
# has_fast_reactions(model)
has_fast_reactions <- function(model) {

  if (!is_azrmod(model))
    stop("has_fast_reactions: input argument is not an AZRmodel")

  if (getNumberOfReactionsAZRmodel(model) < 1) return(FALSE)

  for (k in 1:getNumberOfReactionsAZRmodel(model))
    if (model$reactions[[k]]$fast)
      return(TRUE)

  return(FALSE)
}

###############################################################################
# has_algebraic
###############################################################################
# Checks if the model contains algebraic states
#
# @param model An AZRmodel object
# @return TRUE if model contains algebraic states, FALSE otherwise
# @examples
# model <- exNovakTyson
# has_algebraic(model)
has_algebraic <- function(model) {

  if (!is_azrmod(model))
    stop("has_algebraic: input argument is not an AZRmodel")

  if (getNumberOfAlgebraicAZRmodel(model) < 1) return(FALSE)

  return(TRUE)
}

###############################################################################
# has_constraints
###############################################################################
# Checks if the model contains state constraints
#
# @param model An AZRmodel object
# @return TRUE if model contains state constraints, FALSE otherwise
# @examples
has_constraints <- function(model) {

  if (!is_azrmod(model))
    stop("has_constraints: input argument is not an AZRmodel")

  if (!is.null(attr(model,"originalModel")))
    model <- attr(model,"originalModel")

  for (k in 1:getNumberOfStatesAZRmodel(model)) {
    if (!is.null(model$states[[k]]$lowConstraint))
      return(TRUE)
    if (!is.null(model$states[[k]]$highConstraint))
      return(TRUE)
  }

  return(FALSE)
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
# add_state
# get_state
# set_state
# delete_state
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
# model <- add_state(model,'Cyclin',IC=0.3,'Re1-Re2')
# model <- add_state(model,'Parasites',IC=1e9,'(GR-KR)*Paasites')
# @export
add_state <- function(model,
                             name = NULL,
                             IC = NULL,
                             ODE = NULL,
                             lowConstraint = NULL,
                             highConstraint = NULL,
                             type = NULL,
                             compartment = NULL,
                             unittype = NULL,
                             notes = NULL) {

  if (!is_azrmod(model))
    stop("add_state: model argument is not an AZRmodel")

  if (is.null(name))
    stop("add_state: name is a required input argument")

  if (is.null(IC))
    stop("add_state: IC is a required input argument")

  if (is.null(ODE))
    stop("add_state: ODE is a required input argument")

  if (is.null(lowConstraint) && !is.null(highConstraint))
    stop("add_state: If highConstraint is defined also lowConstraint needs to be defined")

  if (!is.null(lowConstraint) && is.null(highConstraint))
    stop("add_state: If lowConstraint is defined also highConstraint needs to be defined")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("add_state: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("add_state: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

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
# get_state(model,1)
# @export
get_state <- function(model, index) {

  if (!is_azrmod(model))
    stop("get_state: model argument is not an AZRmodel")

  if (getNumberOfStatesAZRmodel(model) < index)
    stop("get_state: value of index larger than number of states.")

  if (index < 1)
    stop("get_state: value of index should be larger than 0.")

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
# model <- set_state(model,1,IC=0.5)
# @export
set_state <- function(model, index,
                             name = NULL,
                             IC = NULL,
                             ODE = NULL,
                             lowConstraint = NULL,
                             highConstraint = NULL,
                             type = NULL,
                             compartment = NULL,
                             unittype = NULL,
                             notes = NULL) {

  if (!is_azrmod(model))
    stop("set_state: model argument is not an AZRmodel")

  if (getNumberOfStatesAZRmodel(model) < index)
    stop("set_state: value of index larger than number of states.")

  if (index < 1)
    stop("set_state: value of index should be larger than 0.")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("set_state: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("set_state: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

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
# delete_state(model,1)
# @export
delete_state <- function(model, index) {

  if (!is_azrmod(model))
    stop("delete_state: model argument is not an AZRmodel")

  if (getNumberOfStatesAZRmodel(model) < index)
    stop("delete_state: value of index larger than number of states.")

  if (index < 1)
    stop("delete_state: value of index should be larger than 0.")

  # check if state has inputs
  if (getNumberOfInputsAZRmodel(model) > 0) {
    inputstates = c()
    for (k in 1:getNumberOfInputsAZRmodel(model)) inputstates = c(inputstates,model$inputs[[k]]$stateindex)
    if (index %in% inputstates)
      stop("delete_state: The state contains inputs. Please delete them first and then the state.")
  }

  model$states[[index]] <- NULL

  return(model)
}


###############################################################################
# Parameter handling functions
# add_parameter
# get_parameter
# set_parameter
# delete_parameter
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
# add_parameter(model,'ka',2,notes="Hello World")
# add_parameter(model,'F',value=0.5)
# @export
add_parameter <- function(model,
                             name = NULL,
                             value = NULL,
                             type = NULL,
                             compartment = NULL,
                             unittype = NULL,
                             notes = NULL,
                             estimate = FALSE,
                             regressor = FALSE) {

  if (!is_azrmod(model))
    stop("add_parameter: model argument is not an AZRmodel")

  if (is.null(name))
    stop("add_parameter: name is a required input argument")

  if (is.null(value))
    stop("add_parameter: value is a required input argument")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("add_parameter: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("add_parameter: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

  if (estimate && regressor)
    stop("add_state: estimate and regressor can not be TRUE at the same time")

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
# get_parameter(model,1)
# @export
get_parameter <- function(model, index) {

  if (!is_azrmod(model))
    stop("get_parameter: model argument is not an AZRmodel")

  if (getNumberOfParametersAZRmodel(model) < index)
    stop("get_parameter: value of index larger than number of states.")

  if (index < 1)
    stop("get_parameter: value of index should be larger than 0.")

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
# set_parameter(model,1,value=0.14)
# @export
set_parameter <- function(model, index,
                                 name = NULL,
                                 value = NULL,
                                 type = NULL,
                                 compartment = NULL,
                                 unittype = NULL,
                                 notes = NULL,
                                 estimate = FALSE,
                                 regressor = FALSE) {

  if (!is_azrmod(model))
    stop("set_parameter: model argument is not an AZRmodel")

  if (getNumberOfParametersAZRmodel(model) < index)
    stop("set_parameter: value of index larger than number of parameters.")

  if (index < 1)
    stop("set_parameter: value of index should be larger than 0.")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("set_parameter: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("set_parameter: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

  if (!is.null(estimate) && !is.null(regressor) && estimate && regressor)
    stop("add_state: estimate and regressor can not be TRUE at the same time")

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
# delete_parameter(model,1)
# @export
delete_parameter <- function(model, index) {

  if (!is_azrmod(model))
    stop("delete_parameter: model argument is not an AZRmodel")

  if (getNumberOfParametersAZRmodel(model) < index)
    stop("delete_parameter: value of index larger than number of parameters.")

  if (index < 1)
    stop("delete_parameter: value of index should be larger than 0.")

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
# add_variable
# get_variable
# set_variable
# delete_variable
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
# add_variable(model,'abc','a+b',notes="hello")
# add_variable(model,'F',formula='a+b')
# @export
add_variable <- function(model,
                                 name = NULL,
                                 formula = NULL,
                                 type = NULL,
                                 compartment = NULL,
                                 unittype = NULL,
                                 notes = NULL) {

  if (!is_azrmod(model))
    stop("add_variable: model argument is not an AZRmodel")

  if (is.null(name))
    stop("add_variable: name is a required input argument")

  if (is.null(formula))
    stop("add_variable: formula is a required input argument")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("add_variable: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("add_variable: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

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
# get_variable(model,1)
# @export
get_variable <- function(model, index) {

  if (!is_azrmod(model))
    stop("get_variable: model argument is not an AZRmodel")

  if (getNumberOfVariablesAZRmodel(model) < index)
    stop("get_variable: value of index larger than number of variables.")

  if (index < 1)
    stop("get_variable: value of index should be larger than 0.")

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
# model <- set_variable(model,1,formula='Cyclin')
# @export
set_variable <- function(model, index,
                                 name = NULL,
                                 formula = NULL,
                                 type = NULL,
                                 compartment = NULL,
                                 unittype = NULL,
                                 notes = NULL) {

  if (!is_azrmod(model))
    stop("set_variable: model argument is not an AZRmodel")

  if (getNumberOfVariablesAZRmodel(model) < index)
    stop("set_variable: value of index larger than number of variables.")

  if (index < 1)
    stop("set_variable: value of index should be larger than 0.")

  if (!is.null(type) && !(type %in% c("isSpecie", "isCompartment", "isParameter")))
    stop("set_variable: wrong definition of 'type'. Needs to be 'isSpecie', 'isCompartment', or 'isParameter'")

  if (!is.null(unittype) && !(unittype %in% c("amount", "concentration")))
    stop("set_variable: wrong definition of 'unittype'. Needs to be 'amount' or 'concentration'")

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
# delete_variable(model,1)
# @export
delete_variable <- function(model, index) {

  if (!is_azrmod(model))
    stop("delete_variable: model argument is not an AZRmodel")

  if (getNumberOfVariablesAZRmodel(model) < index)
    stop("delete_variable: value of index larger than number of variables.")

  if (index < 1)
    stop("delete_variable: value of index should be larger than 0.")

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
# add_reaction
# get_reaction
# set_reaction
# delete_reaction
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
# add_reaction(model,'R1','hello',notes='hello notes',reversible=TRUE)
# add_reaction(model,'R2',formula='a+b')
# @export
add_reaction <- function(model,
                                name = NULL,
                                formula = NULL,
                                notes = NULL,
                                reversible = FALSE,
                                fast = FALSE) {

  if (!is_azrmod(model))
    stop("add_reaction: model argument is not an AZRmodel")

  if (is.null(name))
    stop("add_reaction: name is a required input argument")

  if (is.null(formula))
    stop("add_reaction: formula is a required input argument")

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
# get_reaction(model,19)
# @export
get_reaction <- function(model, index) {

  if (!is_azrmod(model))
    stop("get_reaction: model argument is not an AZRmodel")

  if (getNumberOfReactionsAZRmodel(model) < index)
    stop("get_reaction: value of index larger than number of reactions.")

  if (index < 1)
    stop("get_reaction: value of index should be larger than 0.")

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
# model <- set_reaction(model,19,formula='R1')
# @export
set_reaction <- function(model, index,
                                name = NULL,
                                formula = NULL,
                                notes = NULL,
                                reversible = NULL,
                                fast = NULL) {

  if (!is_azrmod(model))
    stop("set_reaction: model argument is not an AZRmodel")

  if (getNumberOfReactionsAZRmodel(model) < index)
    stop("set_reaction: value of index larger than number of reactions.")

  if (index < 1)
    stop("set_reaction: value of index should be larger than 0.")

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
# delete_reaction(model,1)
# @export
delete_reaction <- function(model, index) {

  if (!is_azrmod(model))
    stop("delete_reaction: model argument is not an AZRmodel")

  if (getNumberOfReactionsAZRmodel(model) < index)
    stop("delete_reaction: value of index larger than number of reactions.")

  if (index < 1)
    stop("delete_reaction: value of index should be larger than 0.")

  model$reactions[[index]] <- NULL

  return(model)
}


###############################################################################
# Function handling functions
# add_function
# get_function
# set_function
# del_function
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
# add_function(model,'ADD',arguments="x,y",formula="x+y")
# add_function(model,'MM',arguments='X,VMAX,KM',formula='VMAX*X/(X+KM)')
# @export
add_function <- function(model,
                                name = NULL,
                                arguments = NULL,
                                formula = NULL,
                                notes = NULL) {

  if (!is_azrmod(model))
    stop("add_function: model argument is not an AZRmodel")

  if (is.null(name))
    stop("add_function: name is a required input argument")

  if (is.null(arguments))
    stop("add_function: arguments is a required input argument")

  if (is.null(formula))
    stop("add_function: formula is a required input argument")

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
# model <- add_function(model,'ADD',arguments="x,y",formula="x+y")
# get_function(model,1)
# @export
get_function <- function(model, index) {

  if (!is_azrmod(model))
    stop("get_function: model argument is not an AZRmodel")

  if (getNumberOfFunctionsAZRmodel(model) < index)
    stop("get_function: value of index larger than number of functions.")

  if (index < 1)
    stop("get_function: value of index should be larger than 0.")

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
# model <- add_function(model,'ADD',arguments="x,y",formula="x+y")
# set_function(model,1,formula='x*y')
# @export
set_function <- function(model, index,
                                name = NULL,
                                arguments = NULL,
                                formula = NULL,
                                notes = NULL) {

  if (!is_azrmod(model))
    stop("set_function: model argument is not an AZRmodel")

  if (getNumberOfFunctionsAZRmodel(model) < index)
    stop("set_function: value of index larger than number of functions.")

  if (index < 1)
    stop("set_function: value of index should be larger than 0.")

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
# model <- add_function(model,'ADD',arguments="x,y",formula="x+y")
# del_function(model,1)
# @export
del_function <- function(model, index) {

  if (!is_azrmod(model))
    stop("del_function: model argument is not an AZRmodel")

  if (getNumberOfFunctionsAZRmodel(model) < index)
    stop("del_function: value of index larger than number of functions.")

  if (index < 1)
    stop("del_function: value of index should be larger than 0.")

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

  if (!is_azrmod(model))
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

  if (!is_azrmod(model))
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

  if (!is_azrmod(model))
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

  if (!is_azrmod(model))
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

  if (!is_azrmod(model))
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

  if (!is_azrmod(model))
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

  if (!is_azrmod(model))
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

  if (!is_azrmod(model))
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

  if (!is_azrmod(model))
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

  if (!is_azrmod(model))
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

  if (!is_azrmod(model))
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

  if (!is_azrmod(model))
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
# model <- add_state(model,'Cyclin',IC=0.3,'Re1-Re2')
# model <- add_state(model,'Cyclin2',IC=0.3,'Re2-Re3')
# model <- addInputAZRmodel(model,stateindex=2)
# model <- addInputAZRmodel(model,factors='F',stateindex=2)
# model <- addInputAZRmodel(model,factors=c('+F','(1-F)'),stateindex=c(1,2))
# @export
addInputAZRmodel <- function(model,
                             stateindex = NULL,
                             factors = "+1") {

  if (!is_azrmod(model))
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
    factors[k] <- strtrimM(factors[k])
    if (substr(factors[k],1,1)!="+") factors[k] <- paste("+",factors[k],sep="")
  }

  terms <- factors
  for (k in 1:length(terms)) terms[k] <- paste(strremWhite(factors[k]),"*",name, sep="")
  model <- add_parameter(model,name=name,value=0)
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
# model <- add_state(model,'Cyclin',IC=0.3,'Re1-Re2')
# model <- add_state(model,'Cyclin2',IC=0.3,'Re2-Re3')
# model <- addInputAZRmodel(model,stateindex=2)
# getInputAZRmodel(model,1)
# @export
getInputAZRmodel <- function(model, index) {

  if (!is_azrmod(model))
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
# model <- add_state(model,'Cyclin',IC=0.3,'Re1-Re2')
# model <- add_state(model,'Cyclin2',IC=0.3,'Re2-Re3')
# model <- addInputAZRmodel(model,stateindex=2)
# delInputAZRmodel(model,1)
# @export
delInputAZRmodel <- function(model, index) {

  if (!is_azrmod(model))
    stop("delInputAZRmodel: model argument is not an AZRmodel")

  if (getNumberOfInputsAZRmodel(model) < index)
    stop("delInputAZRmodel: value of index larger than number of inputs.")

  if (index < 1)
    stop("delInputAZRmodel: value of index should be larger than 0.")

  # Remove INPUT term in ODE
  si <- model$inputs[[index]]$stateindex
  for (k in 1:length(si))
    model$states[[si[k]]]$ODE <- strrepM(strremWhite(model$states[[si[k]]]$ODE),model$inputs[[index]]$terms[[k]],"")

  # Save parindex of input
  parindex <- model$inputs[[index]]$parindex

  # Remove input
  model$inputs[[index]] <- NULL

  # Remove INPUT parameter
  model <- delete_parameter(model,parindex)

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

  if (!is_azrmod(model))
    stop("addOutputAZRmodel: model argument is not an AZRmodel")

  if (is.null(formula))
    stop("addOutputAZRmodel: formula is a required input argument")

  name <- paste("OUTPUT", getNumberOfOutputsAZRmodel(model)+1, sep="")

  model <- add_variable(model,name,formula,notes=notes)

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

  if (!is_azrmod(model))
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

  if (!is_azrmod(model))
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
  model <- delete_variable(model,varindex)

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

  if (!is_azrmod(model))
    stop("getAllStatesAZRmodel: model argument is not an AZRmodel")

  statenames <- c()
  stateICs <- c()
  stateODEs <- c()
  if (getNumberOfStatesAZRmodel(model) > 0) {
    for (k in 1:getNumberOfStatesAZRmodel(model)) {
      x <- get_state(model,k)
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

  if (!is_azrmod(model))
    stop("getAllParametersAZRmodel: model argument is not an AZRmodel")

  paramnames <- c()
  paramvalues <- c()
  paramestimate <- c()
  paramregressor <- c()
  if (getNumberOfParametersAZRmodel(model) > 0) {
    for (k in 1:getNumberOfParametersAZRmodel(model)) {
      x <- get_parameter(model,k)
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

  if (!is_azrmod(model))
    stop("getAllVariablesAZRmodel: model argument is not an AZRmodel")

  varnames <- c()
  varformulas <- c()
  if (getNumberOfVariablesAZRmodel(model) > 0) {
    for (k in 1:getNumberOfVariablesAZRmodel(model)) {
      x <- get_variable(model,k)
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

  if (!is_azrmod(model))
    stop("getAllReactionsAZRmodel: model argument is not an AZRmodel")

  reacnames <- c()
  reacformulas <- c()
  reacfast <- c()
  reacreversible <- c()
  if (getNumberOfReactionsAZRmodel(model) > 0) {
    for (k in 1:getNumberOfReactionsAZRmodel(model)) {
      x <- get_reaction(model,k)
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

# Get information about all functions in an AZRmodel
#
# @param model An AZRmodel
# @return A list with the following components:
# \item{funcnames}{Vector with function names}
# \item{funcformulas}{Vector with function formulas}
# \item{funcarguments}{Flag for fast reactions}
# @examples
# @export
getAllFunctionsAZRmodel <- function(model) {

  if (!is_azrmod(model))
    stop("getAllFunctionsAZRmodel: model argument is not an AZRmodel")

  funcnames <- c()
  funcformulas <- c()
  funcarguments <- c()
  if (getNumberOfFunctionsAZRmodel(model) > 0) {
    for (k in 1:getNumberOfFunctionsAZRmodel(model)) {
      x <- get_function(model,k)
      funcnames[k] <- x$name
      funcformulas[k] <- x$formula
      funcarguments[k] <- x$arguments
    }
  }
  names(funcnames) <- funcnames
  names(funcformulas) <- funcnames
  names(funcarguments) <- funcnames
  return(list(funcnames=funcnames,funcformulas=funcformulas,funcarguments=funcarguments))
}

# Get information about all events in an AZRmodel
#
# @param model An AZRmodel
# @return A list with the following components:
# \item{evenames}{Vector with event names}
# \item{evetriggers}{Vector with event triggers}
# \item{evevariables}{Vector with states/parameters affected by the event assignments}
# \item{eveformulas}{Vector with event assignment expressions}
# @examples
# @export
getAllEventsAZRmodel <- function(model) {

  if (!is_azrmod(model))
    stop("getAllEventsAZRmodel: model argument is not an AZRmodel")

  evenames <- c()
  evetriggers <- c()
  evevariables <- c()
  eveformulas <- c()
  if (getNumberOfEventsAZRmodel(model) > 0) {
    for (k in 1:getNumberOfEventsAZRmodel(model)) {
      x <- getEventAZRmodel(model,k)
      evenames[k] <- x$name
      evetriggers[k] <- x$trigger
      # unlist assignments
      y <- unname(unlist(x$assignment))
      evevariables[k] <- list(y[seq(1,length(y),2)])
      eveformulas[k] <- list(y[seq(2,length(y),2)])
    }
  }
  names(evenames) <- evenames
  names(evetriggers) <- evenames
  names(evevariables) <- evenames
  names(eveformulas) <- evenames
  return(list(evenames=evenames,evetriggers=evetriggers,evevariables=evevariables,eveformulas=eveformulas))
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
  content <- readr::read_file(tempfilename)
  # Exchange strings
  searchStrings <- paste("\\b",origStrings,"\\b",sep="")
  for (k in 1:length(searchStrings)) {
    content <- gsub(searchStrings[k], newStrings[k], content)
  }
  # Save modified model
  readr::write_file(content,tempfilename)
  # Load model (without generation of simulation functions)
  model <- importTxtAZRmodel(create_model(),tempfilename)
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
# Useful to remove "" strings before export of a model.
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
  content <- readr::read_file(tempfilename)
  # Exchange string
  content <- gsub(origString, newString, content, fixed=TRUE)
  # Save modified model
  readr::write_file(content,tempfilename)
  # Load model (without generation of simulation functions)
  model <- importTxtAZRmodel(create_model(),tempfilename)
  # Delete temp file
  unlink(tempfilename)
  # Return model
  return(model)
}
