###############################################################################
###############################################################################
# This file contains auxiliary functions for the import functions from TXT
# and TXTBC files to AZRmodels - specifically handling INPUT and OUTPUT
# definitions
###############################################################################
###############################################################################


###############################################################################
# getOutputs - parse and check outputs in the model
# OUTPUT*:
#       - "*": 1,2,3,4,5, ... sequential indices, starting from 1.
#              no number is allowed to be excluded.
#       - output definitions are only allowed to appear in the model
#         variables section
#       - output RHS expressions should be a single component (due to MONOLIX)
#         state or variable
# model.outputs.name:      output name
# model.outputs.formula:   output formula
# model.outputs.notes:     output notes
# model.outputs.varindex:  index of output in model variables
###############################################################################
getOutputs <- function(model) {

  # Get component information
  states     <- get_all_states(model)
  parameters <- get_all_parameters(model)
  variables  <- get_all_variables(model)
  reactions  <- get_all_reactions(model)

  # Check OUTPUT presence in component names other than variables and error in this case
  if (length(grep("OUTPUT",states$statenames))>0)
    stop("getOutputs: 'OUTPUT' present in a state name. This is not allowed")
  if (length(grep("OUTPUT",parameters$paramnames))>0)
    stop("getOutputs: 'OUTPUT' present in a parameter name. This is not allowed")
  if (length(grep("OUTPUT",reactions$reacnames))>0)
    stop("getOutputs: 'OUTPUT' present in a reaction name. This is not allowed")

  # Check if outputs present in variables
  Noutputs <- length(grep("OUTPUT",variables$varnames))
  if (Noutputs==0) {
    # Not outputs present, return
    return(model)
  }

  # Outputs are present - handle them
  # Assumption is that they are in sequence from 1 to Noutputs
  # Otherwise error!
  for (k in 1:Noutputs) {
    # Try to find OUTPUT"k"
    outputName <- paste("OUTPUT",k,sep="")
    outputVarIndex <- grep(outputName,variables$varnames)
    # Check if found (otherwise error)
    if (length(outputVarIndex)==0)
      stop("getOutputs: Please check OUTPUT* definition in the model. * needs to start from 1 and increment by one for each additional output")
    # Output found
    outputFormula <- model$variables[[outputVarIndex]]$formula
    outputNotes <- model$variables[[outputVarIndex]]$notes
    # Add information to the model
    outputInfo <- list(name=outputName, formula=outputFormula, notes=outputNotes, varindex=outputVarIndex)
    model$outputs[[len_outputs(model)+1]] <- outputInfo
    # Check if output formula is a state or a variable and if not then produce a warning
    test <- cbind(unname(which(outputFormula==states$statenames)),unname(which(outputFormula==variables$varnames)))
    if (length(test)==0)
      warning("getOutput: It is better to assign a single state or a single variable
           to the RHS of an OUTPUT instead of an expression! Important for MONOLIX.")
  }
  return(model)
}


###############################################################################
# getInputs - parse and check inputs in the model
# INPUT*:
#       - "*": 1,2,3,4,5, ...
#       - input definitions are only allowed in differential equations
#         but also input definitions in reactions are allowed. In this case
#         the reaction name in the ODE is replaced by the reaction
#         expression to have the input definition in the ODE. Also, the reaction
#         expression in this case is ONLY ALLOWED TO CONTAIN THE INPUT* TERM
#         possibly with pre-factors - multiplicative! But no parentheses!
#       - a prefactor is allowed, e.g. +1.5*INPUT1 or (k1+k2)*INPUT2
#         the prefactors can be arbitrarily chosen. No postfactors are
#         allowed. Additionally the INPUT* element needs to be a
#         multiplicative factor.
#       - Input terms (INPUT* identifier and prefactor) are only allowed to
#         be added (+) to the differential equation.
#       - If INPUT* is already also present as a parameter, it will be used
#         as defined. If the model does not contain an INPUT* parameter, it
#         is added and set by default to "0". INPUT* is only allowed to be
#         defined as a parameter. Not as a variable and not as a reaction.
#       - Only parameters (and numerical values) are allowed to be used in
#         the definition of INPUT* prefactors. But no states, variables and
#         reactions.
# model.inputs.name:       input name
# model.inputs.factors:    cell-array with input factors
# model.inputs.terms:      cell-array with complete input string (for
#                          simpler removing)
# model.inputs.stateindex: vector with stateindices to which the
#                          input is applied
# model.inputs.parindex:   index of the INPUT* parameter definition in
#                          the IQMmodel (used to remove it when
#                          parameters are written to (e.g.) an MLXTRAN
#                          file).
###############################################################################
getInputs <- function(model) {

  # Handle inputs in reactions
  model <- handleINPUTreactions(model)

  # Check inputs only on states
  if (!inputsOnlyOnStates(model))
    stop("getInputs: INPUT definitions not only on states")

  # Get number of inputs and check if ordering OK
  NINPUTS <- getCheckNumberInputs(model)

  if (NINPUTS==0) return(model)

  # Cycle through inputs and generate the input information
  states <- get_all_states(model)
  for (k in 1:NINPUTS) {
    inputName <- paste("INPUT",k,sep="")
    # Find all state indices for this inputs
    inputStateindex <- grep(paste("\\b",inputName,"\\b",sep=""),states$stateODEs)
    # Find parameter index and add a parameter if not yet present
    parameters <- get_all_parameters(model)
    inputParindex <- grep(paste("\\b",inputName,"\\b",sep=""),parameters$paramnames)
    if (length(inputParindex)==0) {
      # Parameter not yet present - add it
      model <- add_parameter(model,name=inputName,value=0)
      inputParindex <- len_parameters(model)
    }

    # Cycle through stateindex and get the input terms
    inputFactors <- c()
    inputTerms <- c()

    for (k2 in inputStateindex) {
      ODE <- strremWhite(states$stateODEs[[k2]])
      results <- getFactorsTermsInput(model,inputName,ODE)
      inputFactors <- cbind(inputFactors,results$inputFactor)
      inputTerms <- cbind(inputTerms,results$inputTerm)
    }

    # Add information to the model - can not use the addInput function ...
    # need to do "manually"

    inputInfo <- list(name=inputName, factors=inputFactors, terms=inputTerms, stateindex=inputStateindex, parindex=inputParindex)
    model$inputs[[len_inputs(model)+1]] <- inputInfo
  }
  return(model)
}


###############################################################################
# getFactorsTermsInput: get factor and input for given input and given ODE
###############################################################################
getFactorsTermsInput <- function(model,inputName,ODE) {
  ix <- strlocateall(ODE,inputName)
  if (length(ix$start)>1)
    stop("getInputs: Same INPUT more than once on an ODE - not allowed!")
  # get ODE text before and after input
  if (ix$start==1) {
    ODEpre <- ""
  } else {
    ODEpre <- strtrimM(substr(ODE,1,ix$start-1))
  }
  ODEpost = strtrimM(substr(ODE,ix$end+1,nchar(ODE)))

  # 1) the INPUT* identifier needs to be the last element in the input term.
  # => ODEpost needs to start by "+" or "-"
  if (nchar(ODEpost) > 0) {
    if (substr(ODEpost,1,1) != "+" && substr(ODEpost,1,1) != "-") {
      stop("The INPUT* identifier must be the last element in the input term.")
    }
  }
  # and contain as many opening
  # parentheses than closing ones. Because otherwise INPUT* would be in at
  # least one parenthesis. => Not allowed.
  npo <- length(strlocateall(ODEpost,'(')$start)
  npc <- length(strlocateall(ODEpost,')')$start)
  if (npo != npc)
    stop("The INPUT* identifier is not allowed to be inside a parentheses.")

  # 2) The last character in the ODEpre string needs to be a '+' or a '*'. So
  # that in the latter case the part in front of the INPUT* identifier can be
  # understood as a factor. Can also be a "-"
  if (nchar(ODEpre) > 0) {
    if (substr(ODEpre,nchar(ODEpre),nchar(ODEpre)) != '+' &&
        substr(ODEpre,nchar(ODEpre),nchar(ODEpre)) != '-' &&
        substr(ODEpre,nchar(ODEpre),nchar(ODEpre)) != '*') {
      stop("The INPUT term can have a multiplicative pre-term and the whole term needs to be additive!")
    }
  }

  # GET THE FACTORS and TERMS
  inputFactor <- NULL
  inputTerm <- NULL

  # 3) If ODEpre is empty or ODEpre(end) = '+'/'-' then the factor is simply: 1/-1
  if (nchar(ODEpre)==0) {
    inputFactor <- "1"
    inputTerm <- inputName
    return(list(inputFactor=inputFactor, inputTerm=inputTerm))
  }
  if (substr(ODEpre,nchar(ODEpre),nchar(ODEpre)) == "+") {
    inputFactor <- "1"
    inputTerm <- paste("+",inputName,sep="")
    return(list(inputFactor=inputFactor, inputTerm=inputTerm))
  }
  if (substr(ODEpre,nchar(ODEpre),nchar(ODEpre)) == "-") {
    inputFactor <- "-1"
    inputTerm <- paste("-",inputName,sep="")
    return(list(inputFactor=inputFactor, inputTerm=inputTerm))
  }

  # 4) Parse the multiplicative input factor: We search for a '+' or '-'
  # outside of parentheses from right to left in ODEpre. If '-' then error.
  po <- 0
  for (k in nchar(ODEpre):1) {
    if (substr(ODEpre,k,k) == '(') po <- po+1
    if (substr(ODEpre,k,k) == ')') po <- po-1
    if (substr(ODEpre,k,k) == "+" && po==0) break
    if (substr(ODEpre,k,k) == "-" && po==0) break
  }

  # 5) extract the factor
  inputFactor <- substr(ODEpre,k,nchar(ODEpre)-1)
  inputTerm <- paste(inputFactor,"*",inputName,sep="")

  # Check if inputFactor contains state, variable, or reaction names
  stateNames     <- get_all_states(model)$statenames
  variableNames  <- get_all_variables(model)$varnames
  reactionNames  <- get_all_reactions(model)$reacnames

  testNames <- c(stateNames,variableNames,reactionNames)

  for (k in 1:length(testNames)) {
    test <- grep(paste("\\b",testNames[k],"\\b",sep=""),inputFactor)
    if (length(test)!=0)
      stop("getInputs: pre-factors of INPUT* terms are not allowed to depend on states, variables, and reactions!")
  }

  # Return (ITS DONE :))
  return(list(inputFactor=inputFactor, inputTerm=inputTerm))
}


###############################################################################
# getCheckNumberInputs: check if inputs only on states
###############################################################################
getCheckNumberInputs <- function(model) {
  states <- get_all_states(model)
  m      <- gregexpr("INPUT[0-9]+",states$stateODEs,perl=TRUE)
  y      <- regmatches(states$stateODEs,m)
  y      <- unique(unlist(y))
  y      <- sort(as.numeric(strrepM(y,"INPUT","")))

  if (length(y) == 0) return(0)

  # Check if numbers ok
  if (min(y) != 1)
    stop("getInputs: Numbering of INPUT definitions not correct. Has to start with 'INPUT1', continue with 'INPUT2', etc.!")
  if (max(y) != length(y))
    stop("getInputs: Numbering of INPUT definitions not correct. Has to start with 'INPUT1', continue with 'INPUT2', etc.!")
  NINPUTS <- length(y)
  return(NINPUTS)
}


###############################################################################
# inputsOnlyOnStates: check if inputs only on states
###############################################################################
inputsOnlyOnStates <- function(model) {
  variables  <- get_all_variables(model)
  reactions  <- get_all_reactions(model)

  result = TRUE
  if (length(grep("\\bINPUT[0-9]+\\b",variables$varnames))>0) result = FALSE
  if (length(grep("\\bINPUT[0-9]+\\b",variables$varformulas))>0) result = FALSE
  if (length(grep("\\bINPUT[0-9]+\\b",reactions$reacnames))>0) result = FALSE
  if (length(grep("\\bINPUT[0-9]+\\b",reactions$reacformulas))>0) result = FALSE

  return(result)
}

###############################################################################
# handleINPUTreactions: handles inputs in reactions
###############################################################################
handleINPUTreactions <- function(model) {
  # Get component information
  states     <- get_all_states(model)
  parameters <- get_all_parameters(model)
  variables  <- get_all_variables(model)
  reactions  <- get_all_reactions(model)

  # Handle INPUTS in reactions => remove reactions and add to ODEs
  ixReacInputs <- grep("INPUT",reactions$reacformulas)

  if (length(ixReacInputs)==0) {
    # No INPUTS in reactions present
    return(model)
  }

  # Inputs in reactions present - handle them
  for (k in 1:length(ixReacInputs)) {
    reacIndex <- ixReacInputs[k]
    reacName <- reactions$reacnames[reacIndex]
    reacFormula <- reactions$reacformulas[reacIndex]
    # Check formula ... no + or - allowed outside parentheses
    terms <- strexplodePC(reacFormula,"\\+")
    if (length(terms) > 1)
      stop("getInputs: model contains reaction with INPUT as additive term. This is not allowed!")
    terms <- strexplodePC(reacFormula,"\\-")
    if (length(terms) > 1)
      stop("getInputs: model contains reaction with INPUT as additive term. This is not allowed!")
    # Check formula - no parentheses allowed
    terms <- strlocateall(reacFormula,"(")
    if (!is.null(terms$start))
      stop("getInputs: model contains reaction with INPUT and parentheses in formula. Parentheses not allowed in this case!")
    terms <- strlocateall(reacFormula,")")
    if (!is.null(terms$start))
      stop("getInputs: model contains reaction with INPUT and parentheses in formula. Parentheses not allowed in this case!")
    # Now it is ensured that INPUT appears outside parentheses and not in an sum of terms on
    # the RHS of the reaction formula - exchange in ODEs can now happen!
    # Cycle through states and update ODEs
    for (k2 in 1:len_states(model)) {
      model$states[[k2]]$ODE <- strrepM(model$states[[k2]]$ODE,reacName,reacFormula)
    }
    # Check if the reaction name is used in other reactions
    for (k2 in 1:len_reactions(model)) {
      if (length(grep(paste("\\b",reacName,"\\b",sep=""),reactions$reacformulas[k2]))>0)
        stop("getInputs: A reaction with an INPUT is present on the RHS of another reaction. This is not allowed!")
    }
  }
  # Now we can remove the reactions
  reacNamesDelete <- reactions$reacnames[ixReacInputs]
  for (k in 1:length(reacNamesDelete)) {
    # find index
    reactions  <- get_all_reactions(model)
    ixReacInputs <- grep(paste("\\b",reacNamesDelete[k],"\\b",sep=""),reactions$reacnames)
    model <- del_reaction(model,ixReacInputs)
  }
  return(model)
}
