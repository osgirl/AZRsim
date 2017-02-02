###############################################################################
###############################################################################
# This file contains the import functions from TXT files to AZRmodels
###############################################################################
###############################################################################

###############################################################################
# importTxtAZRmodel: imports an AZRmodel from a .txt file
###############################################################################
# Import of of a .txt file to an AZRmodel
#
# A text file with an ODE based representation of an AZRmodel is imported to
# an AZRmodel object. This is an auxiliary function that is called from
# within the AZRmodel function.
#
# @param filename Full path with filename of the .txt file
# @return The imported AZRmodel
# @examples
# importTxtAZRmodel("model.txt")

importTxtAZRmodel <- function(model,filename) {

  ################################################################
  # Preparation
  ################################################################

  # Read model file row by row (already split by row)
  modelText <- fileread(filename,collapserows=FALSE)

  # Remove empty rows
  modelText <- gsub("^\\s+$", "", modelText)
  modelText <- modelText[modelText!=""]

  # Remove commented rows
  modelText <- gsub("^\\s*%.*$", "", modelText)
  modelText <- modelText[modelText!=""]

  ################################################################
  # Cut model text into different sections
  ################################################################

  # Find the starts of the different view data
  name_pos <- grep("^\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* MODEL NAME", modelText)
  notes_pos <- grep("^\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* MODEL NOTES", modelText)
  states_pos <- grep("^\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* MODEL STATES", modelText)
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
    stop("importTxtAZRmodel: There are Matlab functions in the model, please get a copy of IQM tools (http://www.intiquan.com/iqm-tools/).")

  ################################################################
  # Parse the rest
  ################################################################
  model$name <- model_name[1]
  if (length(model_name) > 1) warning("importTxtAZRmodel: model name defined over more than one line. Only first line will be used.")
  model$notes <- paste(strtrim(model_notes),collapse="\n")

  # Handled for TXT models specifically
  model <- getStatesTxt(model,model_states)
  model <- getReactionsTxt(model,model_reactions)

  # Handled for TXT and TXTBC models in the same way
  model <- getParameters(model,model_parameters)
  model <- getVariables(model,model_variables)
  model <- getFunctions(model,model_functions)
  model <- getEvents(model,model_events)

  # Finally input and output need to be checked and handled
  model <- getOutputs(model)
  model <- getInputs(model)

  # Return imported model
  return(model)
}

###############################################################################
# getStatesTxt
###############################################################################
getStatesTxt <- function(model,model_states) {

  # get the rows for ODEs
  ODEtest <- grep("d/dt\\(", model_states)
  # get the rows for ARs
  ARtest <- grep("0 =", model_states)
  # get the rows for ICs
  ICtest <- grep("\\(0\\)", model_states)

  # check if ode definitions are present
  if (length(ODEtest)==0){
    stop('getStatesTxt: The model does not contain any states')
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
    stateString <- strtrim(model_states[ODEtest[k]])

    # Parse comments / notes
    commentInfo <- checkgetNotes(stateString)
    stateString <- commentInfo$main
    notesk      <- commentInfo$comment

    # Check if constraint information is present on a state. Syntax: {constraints:[min,max]}
    stateConstraints <- c(NULL,NULL)
    infoStartConstraints <- grep("\\{constraints:", stateString)

    if (length(infoStartConstraints) > 0) {
      stateString <- strremWhite(stateString)

      tempStart <- regexpr("\\{constraints:", stateString)
      tempEnd <- regexpr("\\]}", stateString)

      # get first bracket after {constraints}
      constraintsString <- strtrim(substr(stateString,(tempStart[1]+13),(tempEnd[tempStart<tempEnd][1])))
      stateString <- strtrim(paste(substr(stateString,1,tempStart[1]-1), substr(stateString,(tempEnd[tempStart<tempEnd][1]),nchar(stateString)-2), sep = ""))
      tempStart2 <- regexpr("\\[", constraintsString)
      tempEnd2 <- regexpr("\\]", constraintsString)
      constraintsString <- substr(constraintsString,(tempStart2[1]+1),(tempEnd2[1]-1))
      stateConstraints <- strexplode(constraintsString,',')
      if (length(stateConstraints) != 2) {
        stop('getStatesTxt: A state-constraint information seems to be wrongly defined')
      }
      stateConstraints[1] <- as.numeric(stateConstraints[1])
      stateConstraints[2] <- as.numeric(stateConstraints[2])
    }

    # Parse SBML related information
    SBMLinfo     <- checkGetSBMLinfo(stateString,"getStatesTxt","state")
    typek        <- SBMLinfo$type
    compartmentk <- SBMLinfo$compartment
    unittypek    <- SBMLinfo$unittype
    stateString  <- SBMLinfo$textString

    # extract the state name
    temp <- regexpr(")", stateString)
    test <- substr(stateString,6,(temp[1]-1))
    # check if state name given
    if (nchar(test) == 0) {
      stop("getStatesTxt: At least on state name in ODE definition is not given.")
    }
    namek <- strremWhite(test)

    # extract the state ODE
    temp <- regexpr("=", stateString)
    test <- substr(stateString,temp[1]+1,nchar(stateString))
    # check if ODE given
    if (nchar(test) == 0) {
      stop("getStatesTxt: At least one RHS of an ODE is not given.")
    }
    # The test string contains now the ODE
    ODEk <- strtrim(test)

    # Add state in model with default IC
    model <- addStateAZRmodel(model,name=namek,IC=0,ODE=ODEk,lowConstraint=stateConstraints[1],highConstraint=stateConstraints[2],type=typek,
                              compartment=compartmentk,unittype=unittypek,notes=notesk)
  }

  ###################
  # PROCESS ARs
  ###################
  # run through the ARs and process them
  if (length(ARtest) > 0) {
    for (k in 1:length(ARtest)) {
      # get each single AR
      ARk <- strtrim(model_states[ARtest[k]])

      # Parse comments / notes
      commentInfo <- checkgetNotes(ARk)
      ARk         <- commentInfo$main
      notesk      <- commentInfo$comment

      # Parse SBML related information
      SBMLinfo     <- checkGetSBMLinfo(ARk,"getStatesTxt","algebraic state")
      typek        <- SBMLinfo$type
      compartmentk <- SBMLinfo$compartment
      unittypek    <- SBMLinfo$unittype
      ARformulak   <- SBMLinfo$textString

      # split rhs in formula and variable name
      terms <- strexplode(ARformulak,':')
      if (length(terms) != 2) {
        ARformulak <- terms[1]
        ARnamek <- NULL # keep it empty
        ARick <- NULL
      } else {
        ARformulak <- strtrim(terms[1])
        ARnamek <- strtrim(terms[2])
        ARick <- 0 # default setting
      }

      # Remove 0 = in formula
      terms <- strexplode(ARformulak,'=')
      if (length(terms) != 2 || as.numeric(terms[1])!=0)
        stop("getStatesTxt: error in algebraic state definition")
      ARformulak <- strtrim(terms[2])

      # add algebraic state to the model
      model <- addAlgebraicAZRmodel(model,name=ARnamek,IC=ARick,formula=ARformulak,type=typek,
                                    compartment=compartmentk,unittype=unittypek,notes=notesk)
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
      # extract the state name
      temp <- regexpr("\\(0\\)", ICString)
      stateName <- strtrim(substr(ICString,1,temp[1]-1))
      # extract the states' initial condition
      temp <- regexpr("=", ICString)
      stateIC <- strtrim(substr(ICString,temp[1]+1,nchar(ICString)))
      found <- FALSE
      # add state ic into model
      ix <- veclocate(getAllStatesAZRmodel(model)$statenames==stateName)
      if (length(ix) != 0) {
        model <- setStateAZRmodel(model,ix,IC=stateIC)
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
        ix <- veclocate(algebraic_names==stateName)
        if (length(ix) != 0) {
          model <- setAlgebraicAZRmodel(model,ix,IC=stateIC)
          found <- TRUE
        }
      }
      # Check if found
      if (!found)
        stop("getStatesTxt: An initial condition is defined for a state that is not present in the model")
    }
  }
  return(model)
}


###############################################################################
# getReactionsTxt
###############################################################################
getReactionsTxt <- function(model,model_reactions) {

  # run through the reactions and process them
  if (!is.null(model_reactions)) {
    for (k in 1:length(model_reactions)) {
      reactionString <- strtrim(model_reactions[k])

      # Parse comments / notes
      commentInfo    <- checkgetNotes(reactionString)
      reactionString <- commentInfo$main
      notesk         <- commentInfo$comment

      # extract the reaction name
      temp <- regexpr("=", reactionString)
      test <- strtrim(substr(reactionString,1,(temp[1]-1)))
      # check if reaction name given
      if (nchar(test) == 0) {
        stop("getReactionsTxt: At least one reaction name not given.")
      }
      namek <- strremWhite(test)

      # extract the reaction expression
      formulak = strtrim(substr(reactionString,(temp+1),nchar(reactionString)))

      # check if the "{reversible}" identifier is present.
      flagInfo       <- checkGetFlag(formulak,"{reversible}")
      formulak       <- flagInfo$textString
      reversibleFlag <- flagInfo$flagPresent

      # check if the "{fast}" identifier is present.
      flagInfo <- checkGetFlag(formulak,"{fast}")
      formulak <- flagInfo$textString
      fastFlag <- flagInfo$flagPresent

      # check if variable expression given
      if (nchar(formulak) == 0) {
        stop("getReactionsTxt: At least one reaction definition not given.")
      }

      # Add reaction to model
      model <- addReactionAZRmodel(model,name=namek,formula=formulak,notes=notesk,reversible=reversibleFlag,fast=fastFlag)
    }
  }
  return(model)
}
