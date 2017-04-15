###############################################################################
###############################################################################
# This file contains the import functions from TXTBC files to AZRmodels
###############################################################################
###############################################################################

###############################################################################
# importTxtBcAZRmodel: imports an AZRmodel from a .txtbc file
###############################################################################
# Import of of a .txtbc file to an AZRmodel
#
# A text file with an TEXTBC based representation of an AZRmodel is imported to
# an AZRmodel object. This is an auxiliary function that is called from
# within the AZRmodel function.
#
# @param filename Full path with filename of the .txtbc file
# @return The imported AZRmodel
# @examples
# importTxtBcAZRmodel("model.txtbc")

importTxtBcAZRmodel <- function(model,filename) {

  ################################################################
  # Preparation
  ################################################################

  modelText <- readr::read_lines(filename)

  # Remove empty rows
  modelText <- gsub("^\\s+$", "", modelText)
  modelText <- modelText[modelText!=""]

  # Remove commented rows
  modelText <- gsub("^\\s*#.*$", "", modelText)
  modelText <- modelText[modelText!=""]

  ################################################################
  # Cut model text into different sections
  ################################################################

  # Find the starts of the different view data
  name_pos <- grep("^\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* MODEL NAME", modelText)
  notes_pos <- grep("^\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* MODEL NOTES", modelText)
  states_pos <- grep("^\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* MODEL STATE INFORMATION", modelText)
  parameters_pos <- grep("^\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* MODEL PARAMETERS", modelText)
  variables_pos <- grep("^\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* MODEL VARIABLES", modelText)
  reactions_pos <- grep("^\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* MODEL REACTIONS", modelText)
  functions_pos <- grep("^\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* MODEL FUNCTIONS", modelText)
  events_pos <- grep("^\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* MODEL EVENTS", modelText)
  matlab_functions_pos <- grep("^\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* MODEL MATLAB FUNCTIONS", modelText)

  # Cut out the different pieces and assign them to the modelTextStructure structure
  model_name <- modelText[name_pos:(notes_pos-1)]
  model_notes <- modelText[notes_pos:(states_pos-1)]
  model_states <- modelText[states_pos:(parameters_pos-1)]
  model_parameters <- modelText[parameters_pos:(variables_pos-1)]
  model_variables <- modelText[variables_pos:(reactions_pos-1)]
  model_reactions <- modelText[reactions_pos:(functions_pos-1)]
  model_functions <- modelText[functions_pos:(events_pos-1)]
  if (length(matlab_functions_pos) != 0) {
    # MATLAB function separator present
    model_events <- modelText[events_pos:(matlab_functions_pos-1)]
    model_MATLABfunctions <- modelText[matlab_functions_pos:length(modelText)]
  } else {
    # No MATLAB function separator present
    model_events <- modelText[events_pos:length(modelText)]
    model_MATLABfunctions <- "********** MODEL MATLAB FUNCTIONS"
  }

  # Check if information is present
  if (length(model_name) == 1) { model_name <- "Please provide a model name (one line)" } else { model_name <- model_name[2:length(model_name)] }
  if (length(model_notes) == 1) { model_notes <- "Please provide some model information\nMultiple lines possible" } else { model_notes <- model_notes[2:length(model_notes)] }
  if (length(model_states) == 1) { model_states <- NULL } else { model_states <- model_states[2:length(model_states)] }
  if (length(model_parameters) == 1) { model_parameters <- NULL } else { model_parameters <- model_parameters[2:length(model_parameters)] }
  if (length(model_variables) == 1) { model_variables <- NULL } else { model_variables <- model_variables[2:length(model_variables)] }
  if (length(model_reactions) == 1) { model_reactions <- NULL } else { model_reactions <- model_reactions[2:length(model_reactions)] }
  if (length(model_functions) == 1) { model_functions <- NULL } else { model_functions <- model_functions[2:length(model_functions)] }
  if (length(model_events) == 1) { model_events <- NULL } else { model_events <- model_events[2:length(model_events)] }
  if (length(model_MATLABfunctions) == 1) { model_MATLABfunctions <- NULL } else { model_MATLABfunctions <- model_MATLABfunctions[2:length(model_MATLABfunctions)] }

  # Check if content in model_MATLABfunctions
  if (length(model_MATLABfunctions) > 0)
    stop("importTxtBcAZRmodel: There are Matlab functions in the model, please get a copy of IQM tools (http://www.intiquan.com/iqm-tools/).")

  ################################################################
  # Parse the rest
  ################################################################
  model$name <- model_name[1]
  if (length(model_name) > 1) warning("importTxtBcAZRmodel: model name defined over more than one line. Only first line will be used.")
  model$notes <- paste(strtrimM(model_notes),collapse="\n")
  # Parameters, Variables, Functions, and Events are handled exactly like for TXT model import
  model <- getParameters(model,model_parameters)
  model <- getVariables(model,model_variables)
  model <- getFunctions(model,model_functions)
  model <- getEvents(model,model_events)

  # Only states and reactions are handled differently
  model <- getStatesReactionsTxtBc(model,model_states,model_reactions)

  # Finally input and output need to be checked and handled
  model <- getOutputs(model)
  model <- getInputs(model)

  # Return imported model
  return(model)
}

###############################################################################
# getReactionTerms
###############################################################################
getReactionTerms <- function(reactionPart) {
  if (reactionPart=="") {
    result <- list(names=NULL, factors=NULL)
    return(result)
  }
  # explode into species and stoichiometric coefficients
  # expect the terms to be of the format:
  # factor*species + factor*species + ...
  allTerms <- strexplode(reactionPart,"\\+")

  # check the syntax of the single terms (name or numeric*name)
  reactiontermsNames <- c()
  reactiontermsFactors <- c()

  for (k in 1:length(allTerms)) {
    checkTerms <- strexplode(allTerms[k],"\\*")
    if (length(checkTerms) == 1) {
      reactiontermsNames <- c(reactiontermsNames, checkTerms[1])
      reactiontermsFactors <- c(reactiontermsFactors, 1)
    } else {
      if (length(checkTerms) == 2) {
        reactiontermsNames <- c(reactiontermsNames, checkTerms[2])
        reactiontermsFactors <- c(reactiontermsFactors, checkTerms[1])
      } else {
        stop("getReactionTerms: Syntax error in a reaction equation.")
      }
    }
  }
  result <- list(names=reactiontermsNames, factors=reactiontermsFactors)
  return(result)
}

###############################################################################
# getStatesReactionsTxtBc
###############################################################################
getStatesReactionsTxtBc <- function(model,model_states,model_reactions) {

  ###############################################
  # PARSING OF REACTIONS
  ###############################################

  reactionsList <- list()
  allSpecies <- c()
  k <- 1
  nReac <- 0
  while (k < length(model_reactions)) {
    # Get next row
    rowText = model_reactions[k]
    # Check if row contains a new reaction definition
    if (!is.null(strlocateall(rowText,":")$start)) {
      nReac <- nReac+1
      # It is a new reaction definition
      # Get notes and definition
      reacInfo <- checkgetNotes(rowText)
      reacDef <- strremWhite(reacInfo$main)
      reacNotes <- reacInfo$comment
      # Check FAST flag
      flagCheck <- checkGetFlag(reacDef,"{fast}")
      reacFast <- flagCheck$flagPresent
      reacDef <- flagCheck$textString
      # Get name, LHS, RHS, and reversibility information
      if (!is.null(strlocateall(reacDef,"<=>")$start)) {
        reacReversible <- TRUE
        reacDef <- strrepM(reacDef,"<=>","&&&")
      } else {
        if (!is.null(strlocateall(reacDef,"=>")$start)) {
          reacReversible <- FALSE
          reacDef <- strrepM(reacDef,"=>","&&&")
        } else {
          stop("getStatesReactionsTxtBc: wrong reaction expression")
        }
      }
      # Get name and equation
      terms <- strexplode(reacDef,":")
      reacName <- terms[2]
      reacDef <- terms[1]
      # Get LHS and RHS
      terms <- strexplode(paste(" ",reacDef," ",sep=""),"&&&")
      reacLHS <- strtrimM(terms[1])
      reacRHS <- strtrimM(terms[2])

      # Now parse the reaction kinetics
      # Next row should be "vf="
      k <- k+1
      if (k>length(model_reactions)) stop("getStatesReactionsTxtBc: error in reaction expression")
      rowKineticsVF <- strremWhite(model_reactions[k])
      if (is.null(strlocateall(rowKineticsVF,"vf=")$start))
        stop("getStatesReactionsTxtBc: error in reaction expression")
      rowKineticsVF <- strrepM(rowKineticsVF,"vf=","")
      if (reacReversible) {
        # Next row should be "vr="
        k <- k+1
        if (k>length(model_reactions)) stop("getStatesReactionsTxtBc: error in reaction expression")
        rowKineticsVR <- strremWhite(model_reactions[k])
        if (is.null(strlocateall(rowKineticsVR,"vr=")$start))
          stop("getStatesReactionsTxtBc: error in reaction expression")
        rowKineticsVR <- strrepM(rowKineticsVR,"vr=","")
        reacFormula <- paste(rowKineticsVF,"-",rowKineticsVR,sep="")
      } else {
        rowKineticsVR <- NULL
        reacFormula <- rowKineticsVF
      }
    }

    # Generate model reaction information
    model <- addReactionAZRmodel(model,name=reacName,formula=reacFormula,notes=reacNotes,
                                 reversible=reacReversible,fast=reacFast)

    # Get substrate and product information - needed for ODE construction
    substrateInfo <- getReactionTerms(reacLHS)
    productInfo <- getReactionTerms(reacRHS)

    # Collect all species names
    if (!is.null(substrateInfo$names)) allSpecies <- unique(c(allSpecies, substrateInfo$names))
    if (!is.null(productInfo$names)) allSpecies <- unique(c(allSpecies, productInfo$names))


    reacInfo <- list(name=reacName, substrateNames=substrateInfo$names,
                     substrateFactors=substrateInfo$factors,
                     productNames=productInfo$names,
                     productFactors=productInfo$factors)
    reactionsList[[nReac]] <- reacInfo

    # Increment
    k <- k+1
  }

  ###############################################
  # GENERATION OF STATE INFORMATION
  ###############################################

  # - check the list of species against variables and parameters to
  #   determine the states.
  # - check unittypes ........ BIG THING!!!
  #   define the state substructure
  allSpeciesStates <- c()
  allParNames <- getAllParametersAZRmodel(model)$paramnames
  allVarNames <- getAllVariablesAZRmodel(model)$varnames
  if (length(allSpecies) > 0) {
    for (k in 1:length(allSpecies)) {
      parameterIndex <- strmatch(allSpecies[k],allParNames)
      variableIndex <- strmatch(allSpecies[k],allVarNames)
      if (is.null(parameterIndex) && is.null(variableIndex))
        allSpeciesStates <- c(allSpeciesStates, allSpecies[k])
    }
  }

  # Now generate ODE information for each species defined by reaction expresssions and
  # use standard parsing of states from TXT model
  for (k in 1:length(allSpeciesStates)) {
    stateName <- allSpeciesStates[k]
    # Construct RHS of ODE
    stateODE <- ""
    for (k2 in 1:length(reactionsList)) {
      reacName <- reactionsList[[k2]]$name
      # Check substrate names
      ix <- strmatch(stateName,reactionsList[[k2]]$substrateNames)
      if (!is.null(ix))
        stateODE <- paste(stateODE,"-",reactionsList[[k2]]$substrateFactors[ix],"*",reacName,sep="")
      # Check product names
      ix <- strmatch(stateName,reactionsList[[k2]]$productNames)
      if (!is.null(ix))
          stateODE <- paste(stateODE,"+",reactionsList[[k2]]$productFactors[ix],"*",reacName,sep="")
    }
    ODEtext <- paste("d/dt(",stateName,") = ",stateODE,sep="")

    # Add to model_states
    model_states <- c(ODEtext,model_states)
  }

  ###############################################
  # PARSING OF STATE INFORMATION
  ###############################################
  model <- getStatesTxtBc(model,model_states)

  # Ready
  return(model)
}

###############################################################################
# getStatesTxtBc
###############################################################################
getStatesTxtBc <- function(model,model_states) {

  # get the rows for ODEs
  ODEtest <- grep("d/dt\\(", model_states)
  # get the rows for ARs
  ARtest <- grep("0 =", model_states)
  # get the rows for ICs
  ICtest <- grep("\\(0\\)", model_states)

  # check if ode definitions are present
  if (length(ODEtest)==0){
    stop('getStatesTxtBc: The model does not contain any states')
  }

  # check if ODEs, ARs, and ICs come subsequently
  if (length(ICtest) != 0) {
    if (max(ODEtest) > min(ICtest)) {
      stop('Initial conditions have to be defined after the definition of the ODEs.')
    }
  }
  if (length(ARtest) != 0) {
    if (max(ODEtest) > min(ARtest)) {
      stop('Algebraic rules have to be defined after the definition of the ODEs.')
    }
  }
  if ((length(ARtest) != 0) & (length(ICtest) != 0)) {
    if (max(ARtest) > min(ICtest)) {
      stop('Initial conditions have to be defined after the definition of the algebraic rules.')
    }
  }

  ###################
  # PROCESS ODEs
  ###################
  for (k in 1:length(ODEtest)) {
    stateString <- strtrimM(model_states[ODEtest[k]])

    # extract the state name
    temp <- regexpr(")", stateString)
    test <- substr(stateString,6,(temp[1]-1))
    # check if state name given
    if (nchar(test) == 0) {
      stop("getStatesTxtBc: At least on state name in ODE definition is not given.")
    }
    namek <- strremWhite(test)

    # extract the state ODE
    temp <- regexpr("=", stateString)
    test <- substr(stateString,temp[1]+1,nchar(stateString))
    # check if ODE given
    if (nchar(test) == 0) {
      stop("getStatesTxtBc: At least one RHS of an ODE is not given.")
    }
    # The test string contains now the ODE
    ODEk <- strtrimM(test)

    # Add state in model with default IC
    model <- addStateAZRmodel(model,name=namek,IC=0,ODE=ODEk)
  }

  ###################
  # PROCESS ARs
  ###################
  # run through the ARs and process them
  if (length(ARtest) > 0) {
    for (k in 1:length(ARtest)) {
      # get each single AR
      ARk <- strtrimM(model_states[ARtest[k]])

      # split rhs in formula and variable name
      terms <- strexplode(ARk,':')
      if (length(terms) != 2) {
        ARformulak <- terms[1]
        ARnamek <- NULL # keep it empty
        ARick <- NULL
      } else {
        ARformulak <- strtrimM(terms[1])
        ARnamek <- strtrimM(terms[2])
        ARick <- 0 # default setting
      }

      # Remove 0 = in formula
      terms <- strexplode(ARformulak,'=')
      if (length(terms) != 2 || as.numeric(terms[1])!=0)
        stop("getStatesTxtBc: error in algebraic state definition")
      ARformulak <- strtrimM(terms[2])

      # add algebraic state to the model
      model <- addAlgebraicAZRmodel(model,name=ARnamek,IC=ARick,formula=ARformulak)
    }
  }

  ###################
  # PROCESS ICs
  ###################
  # run through the initial conditions and add them they can have a different order than the odes. if an initial
  # condition is not defined for a certain state then it is set to zero by default
  # First check if any initial conditions are given - if not then don't execute this part!
  if (length(ICtest) > 0) {
    for (k1 in 1:length(ICtest)) {
      ICString <- strremWhite(model_states[ICtest[k1]])

      # Parse comments / notes
      commentInfo <- checkgetNotes(ICString)
      ICString    <- commentInfo$main
      notesk      <- commentInfo$comment

      # Check if constraint information is present on a state. Syntax: {constraints:[min,max]}
      stateConstraints <- c(NULL,NULL)
      infoStartConstraints <- grep("\\{constraints:", ICString)

      if (length(infoStartConstraints) > 0) {
        ICString <- strremWhite(ICString)

        tempStart <- regexpr("\\{constraints:", ICString)
        tempEnd <- regexpr("\\]\\}", ICString)

        # get first bracket after {constraints}
        constraintsString <- strtrimM(substr(ICString,(tempStart[1]+13),(tempEnd[tempStart<tempEnd][1])))
        ICString <- strtrimM(paste(substr(ICString,1,tempStart[1]-1), substr(ICString,(tempEnd[tempStart<tempEnd][1]),nchar(ICString)-2), sep = ""))
        tempStart2 <- regexpr("\\[", constraintsString)
        tempEnd2 <- regexpr("\\]", constraintsString)
        constraintsString <- substr(constraintsString,(tempStart2[1]+1),(tempEnd2[1]-1))
        stateConstraints <- strexplode(constraintsString,',')
        if (length(stateConstraints) != 2) {
          stop('getStatesTxtBc: A state-constraint information seems to be wrongly defined')
        }
        stateConstraints[1] <- as.numeric(stateConstraints[1])
        stateConstraints[2] <- as.numeric(stateConstraints[2])
      }

      # Parse SBML related information
      SBMLinfo     <- checkGetSBMLinfo(ICString,"getStatesTxtBc","state")
      typek        <- SBMLinfo$type
      compartmentk <- SBMLinfo$compartment
      unittypek    <- SBMLinfo$unittype
      ICString     <- SBMLinfo$textString

      # extract the state name
      temp <- regexpr("\\(0\\)", ICString)
      stateName <- strtrimM(substr(ICString,1,temp[1]-1))

      # extract the states' initial condition
      temp <- regexpr("=", ICString)
      stateIC <- strtrimM(substr(ICString,temp[1]+1,nchar(ICString)))

      # add state information into model
      found <- FALSE

      ix <- unname(which(getAllStatesAZRmodel(model)$statenames==stateName))

      if (length(ix) > 1)
        stop("getStatesTxtBc: error in model definition - a state appears more than once.")

      if (length(ix) != 0) {
        model <- setStateAZRmodel(model,ix,IC=stateIC,lowConstraint=stateConstraints[1],
                                  highConstraint=stateConstraints[2],type=typek,
                                  compartment=compartmentk,unittype=unittypek,notes=notesk)
        found <- TRUE
      }

      # add algebraic ic into model
      if (!found && getNumberOfAlgebraicAZRmodel(model)>0) {
        algebraic_names = c()
        for (k2 in 1:getNumberOfAlgebraicAZRmodel(model)) {
          if (!is.null(model$algebraic[[k2]]$name)) {
            algebraic_names <- cbind(algebraic_names,model$algebraic[[k2]]$name)
          } else {
            algebraic_names <- cbind(algebraic_names,"UNDEFINED_AR_NAME")
          }
        }
        ix <- unname(which(algebraic_names==stateName))
        if (length(ix) != 0) {
          model <- setAlgebraicAZRmodel(model,ix,IC=stateIC,type=typek,compartment=compartmentk,unittype=unittypek,
                                        notes=notesk)
          found <- TRUE
        }
      }
      # Check if found
      if (!found)
        stop("getStatesTxtBc: An initial condition is defined for a state that is not present in the model")
    }
  }
  return(model)
}
