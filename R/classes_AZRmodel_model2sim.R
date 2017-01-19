###############################################################################
###############################################################################
# Generation of deSolve related simulation functions for AZRmodels.
###############################################################################
###############################################################################


###############################################################################
# genSimFunctions: Generate simulation functions for deSolve
#                  and attach them to the model
###############################################################################
# Generate a simulation function for deSolve
#
# Several functions will be generated and attached to the model as attributes
# ODEsim: ODE simulation function for deSolve
# VARsim: function returning variable and reaction values
# EVEsim: Event trigger function
# EASsim: Event assignment function
# Additionally, before generating the simulation functions, the model will be
# updated to be able to handle constraints (only if constraints are present)
#
# In the case that the model contains DAEs, the following function will be generated
# DAEsim: DAE simulation function for deSolve
genSimFunctions <- function (model) {

  if (!is.AZRmodel(model))
    stop("genSimFunctions: input argument is not an AZRmodel")

  if (getNumberOfStatesAZRmodel(model)==0)
    stop("genSimFunctions: model does not contain any dynamic states")

  ##############################################################################
  # Update math in the model based on components
  ##############################################################################

  # Handle potential constraints on states
  model <- handleConstraintsSim(model)

  # Handle input math on models
  model <- implementALLinputMath(model)

  ##############################################################################
  # Handle deSolve specific things in the model
  #   This is not done in the model itself but in a copy that is passed
  #   to the simulation function generation function
  ##############################################################################

  # Replace names of specific functions and kinetic rate laws by adding the AZRsim:::
  # namespace identifier in front.
  origStrings <- c("gt","ge","lt","le","mod","and","or","multiply","piecewise",
                   "interp0","interp1","interpcs","interpcse")
  origStrings <- c(origStrings, getAllKineticRateLaws())
  newStrings <- paste("AZRsim:::",origStrings,sep="")
  modelDeSolveSpecific <- renameElementsAZRmodel(model, origStrings, newStrings)

  # Replace syntax for interpolation functions from
  # interpx([comma sep vector],[comma sep vector],element) to
  # interpx(c(comma sep vector),c(comma sep vector),element)
  # More in general all "vector definitions" with "[...]" are replaced by "c(...)"
  modelDeSolveSpecific <- handleVectorSyntaxDeSolve(modelDeSolveSpecific)


  ##############################################################################
  # Generate all simulation functions for deSolve
  ##############################################################################

  # Generate and add ODEsim function
  attr(model,"ODEsim") <- genODEsim(modelDeSolveSpecific)

  # Generate and add VARsim function
  attr(model,"VARsim") <- genVARsim(modelDeSolveSpecific)

  # Generate and add ROOTsim function if events present
  if (getNumberOfEventsAZRmodel(model) > 0) {
    attr(model,"ROOTsim") <- genROOTsim(modelDeSolveSpecific)
    attr(model,"EVASsim") <- genEVASsim(modelDeSolveSpecific)
  } else {
    attr(model,"ROOTsim") <- NULL
    attr(model,"EVASsim") <- NULL
  }

  # Generate and add non-numerical IC handling function
  attr(model,"nnICsim") <- nnICsim(modelDeSolveSpecific)

  ##############################################################################
  # Return the updated model
  ##############################################################################
  return(model)
}


###############################################################################
# handleVectorSyntaxDeSolve: Handles vector syntax
###############################################################################
handleVectorSyntaxDeSolve <- function(model) {

  # Handle ODEs
  for (k in 1:getNumberOfStatesAZRmodel(model)) {
    formula <- model$states[[k]]$ODE
    # Check if square brackets present in formula - if so then lets
    # make R vectors out of it
    if (!is.null(AZRaux:::strlocateall(formula,"[")$start)) {
      # Square brackets present
      formula <- AZRaux:::strremWhite(formula)
      formula <- AZRaux:::strrep(formula,"[","c(")
      formula <- AZRaux:::strrep(formula,"]",")")
      model$states[[k]]$ODE <- formula
    }
  }

  # Handle Variables
  if (getNumberOfVariablesAZRmodel(model) > 0) {
    for (k in 1:getNumberOfVariablesAZRmodel(model)) {
      formula <- model$variables[[k]]$formula
      # Check if square brackets present in formula - if so then lets
      # make R vectors out of it
      if (!is.null(AZRaux:::strlocateall(formula,"[")$start)) {
        # Square brackets present
        formula <- AZRaux:::strremWhite(formula)
        formula <- AZRaux:::strrep(formula,"[","c(")
        formula <- AZRaux:::strrep(formula,"]",")")
        model$variables[[k]]$formula <- formula
      }
    }
  }

  # Handle Reactions
  if (getNumberOfReactionsAZRmodel(model) > 0) {
    for (k in 1:getNumberOfReactionsAZRmodel(model)) {
      formula <- model$reactions[[k]]$formula
      # Check if square brackets present in formula - if so then lets
      # make R vectors out of it
      if (!is.null(AZRaux:::strlocateall(formula,"[")$start)) {
        # Square brackets present
        formula <- AZRaux:::strremWhite(formula)
        formula <- AZRaux:::strrep(formula,"[","c(")
        formula <- AZRaux:::strrep(formula,"]",")")
        model$reactions[[k]]$formula <- formula
      }
    }
  }

  return(model)
}


###############################################################################
# handleConstraintsSim: Handles constraints in an AZRmodel
###############################################################################
handleConstraintsSim <- function (model) {

  # Cycle through states and check for constraints
  for (k in 1:getNumberOfStatesAZRmodel(model)) {
    si <- getStateAZRmodel(model,k)
    if (!is.null(si$lowConstraint) && !is.null(si$highConstraint)) {
      # Add switch condition variable
      varswitchname <- paste("switchConstraint_",si$name,sep="")
      varswitchcon <- paste("and(ge(",si$name,",",si$lowConstraint,"),le(",si$name,",",si$highConstraint,"))",sep="")
      varswitchformula <- paste("piecewise(1,",varswitchcon,",0)",sep="")
      model <- addVariableAZRmodel(model,name = varswitchname,formula=varswitchformula)
      # Update RHS of ODE and remove constraint information
      ODEswitch <- paste(varswitchname,"*(",si$ODE,")",sep="")
      model <- setStateAZRmodel(model,k,ODE=ODEswitch)
      model$states[[k]]$lowConstraint <- NULL
      model$states[[k]]$highConstraint <- NULL
    }
  }
  # Return the updated model
  return(model)
}

###############################################################################
# nnICsim: Handle default non-numeric initial conditions
###############################################################################
nnICsim <- function (model) {

  nnICfctText <- "function () {\n"

  if (length(model$functions)>0) {
    for (k in 1:length(model$functions)) {
      nnICfctText <- paste(nnICfctText,"  ",model$functions[[k]]$name," <- function(",model$functions[[k]]$arguments,") { ",model$functions[[k]]$formula," }\n",sep="")
    }
    nnICfctText <- paste(nnICfctText,"\n",sep="")
  }

  for (k in 1:length(model$states)) {
    nnICfctText <- paste(nnICfctText,"  try(",model$states[[k]]$name," <- ",model$states[[k]]$IC,", silent=TRUE)\n",sep="")
  }
  nnICfctText <- paste(nnICfctText,"\n",sep="")

  if (length(model$parameters)>0) {
    for (k in 1:length(model$parameters)) {
      nnICfctText <- paste(nnICfctText,"  try(",model$parameters[[k]]$name," <- ",model$parameters[[k]]$value,", silent=TRUE)\n",sep="")
    }
    nnICfctText <- paste(nnICfctText,"\n",sep="")
  }

  if (length(model$variables)>0) {
    for (k in 1:length(model$variables)) {
      nnICfctText <- paste(nnICfctText,"  try(",model$variables[[k]]$name," <- ",model$variables[[k]]$formula,", silent=TRUE)\n",sep="")
    }
    nnICfctText <- paste(nnICfctText,"\n",sep="")
  }

  if (length(model$reactions)>0) {
    for (k in 1:length(model$reactions)) {
      nnICfctText <- paste(nnICfctText,"  try(",model$reactions[[k]]$name," <- ",model$reactions[[k]]$formula,", silent=TRUE)\n",sep="")
    }
    nnICfctText <- paste(nnICfctText,"\n",sep="")
  }

  nnICfctText <- paste(nnICfctText,"  out = c()\n",sep="")

  for (k in 1:length(model$states)) {
    nnICfctText <- paste(nnICfctText,"  out['",model$states[[k]]$name,"'] <- tryCatch(",model$states[[k]]$IC,",error=function(cond) return(NA))\n",sep="")
  }

  nnICfctText <- paste(nnICfctText,"\n",sep="")

  nnICfctText <- paste(nnICfctText,"  if (length(AZRaux:::veclocate(is.na(out)))>0)\n",sep="")
  nnICfctText <- paste(nnICfctText,"    stop('nnICsim: problem in evaluation of non-numerical initial conditions')\n",sep="")

  nnICfctText <- paste(nnICfctText,"\n",sep="")
  nnICfctText <- paste(nnICfctText,"  return(out)\n",sep="")

  nnICfctText <- paste(nnICfctText,"}\n",sep="")

  outFunc <- eval(parse(text=nnICfctText))
}

###############################################################################
# genODEsim: Generate ODEsim function for deSolve
###############################################################################
genODEsim <- function (model) {

  SIMfctText <- "function (time, states, paramvalues) {\n"

  SIMfctText <- paste(SIMfctText,getODEtext(model),sep="")

  SIMfctText <- paste(SIMfctText,"  out <- list(c(\n",sep="")
  if (length(model$states) > 1) {
    for (k in 1:(length(model$states)-1)) {
      SIMfctText <- paste(SIMfctText,"    ddt_",model$states[[k]]$name,",\n",sep="")
    }
  }
  SIMfctText <- paste(SIMfctText,"    ddt_",model$states[[length(model$states)]]$name,"))\n",sep="")
  SIMfctText <- paste(SIMfctText,"\n",sep="")

  SIMfctText <- paste(SIMfctText,"  return(out)\n",sep="")
  SIMfctText <- paste(SIMfctText,"}\n",sep="")

  outFunc <- eval(parse(text=SIMfctText))
}

###############################################################################
# genODEtext: Generate auxiliary text that is often used
###############################################################################
getODEtext <- function(model) {

  ODEtext <- ""

  for (k in 1:length(model$states)) {
    ODEtext <- paste(ODEtext,"  ",model$states[[k]]$name," = states[",k,"]\n",sep="")
  }
  ODEtext <- paste(ODEtext,"\n",sep="")

  if (length(model$parameters)>0) {
    for (k in 1:length(model$parameters)) {
      ODEtext <- paste(ODEtext,"  ",model$parameters[[k]]$name," = paramvalues[",k,"]\n",sep="")
    }
    ODEtext <- paste(ODEtext,"\n",sep="")
  }

  if (length(model$functions)>0) {
    for (k in 1:length(model$functions)) {
      ODEtext <- paste(ODEtext,"  ",model$functions[[k]]$name," <- function(",model$functions[[k]]$arguments,") { ",model$functions[[k]]$formula," }\n",sep="")
    }
    ODEtext <- paste(ODEtext,"\n",sep="")
  }

  if (length(model$variables)>0) {
    for (k in 1:length(model$variables)) {
      ODEtext <- paste(ODEtext,"  ",model$variables[[k]]$name," = ",model$variables[[k]]$formula,"\n",sep="")
    }
    ODEtext <- paste(ODEtext,"\n",sep="")
  }

  if (length(model$reactions)>0) {
    for (k in 1:length(model$reactions)) {
      ODEtext <- paste(ODEtext,"  ",model$reactions[[k]]$name," = ",model$reactions[[k]]$formula,"\n",sep="")
    }
    ODEtext <- paste(ODEtext,"\n",sep="")
  }

  for (k in 1:length(model$states)) {
    ODEtext <- paste(ODEtext,"  ddt_",model$states[[k]]$name," = ",model$states[[k]]$ODE,"\n",sep="")
  }
  ODEtext <- paste(ODEtext,"\n",sep="")
  return(ODEtext)
}


###############################################################################
# genVARsim: Generate VARsim function
###############################################################################
genVARsim <- function (model) {

  VARfctText <- "function (simresODE, paramvalues) {\n"

  VARfctText <- paste(VARfctText,"  simresVAR <- c()\n",sep="")

  if (length(model$variables) > 0 || length(model$reactions) > 0) {
    VARfctText <- paste(VARfctText,"  for (ksimresODE in 1:nrow(simresODE)) {\n",sep="")

    VARfctText <- paste(VARfctText,"    time = simresODE[ksimresODE,1]\n",sep="")

    for (k in 1:length(model$states)) {
      VARfctText <- paste(VARfctText,"    ",model$states[[k]]$name," = simresODE[ksimresODE,",k+1,"]\n",sep="")
    }

    if (length(model$parameters)>0) {
      for (k in 1:length(model$parameters)) {
        VARfctText <- paste(VARfctText,"    ",model$parameters[[k]]$name," = paramvalues[",k,"]\n",sep="")
      }
    }

    if (length(model$functions)>0) {
      for (k in 1:length(model$functions)) {
        VARfctText <- paste(VARfctText,"  ",model$functions[[k]]$name," <- function(",model$functions[[k]]$arguments,") { ",model$functions[[k]]$formula," }\n",sep="")
      }
      VARfctText <- paste(VARfctText,"\n",sep="")
    }

    if (length(model$variables)>0) {
      for (k in 1:length(model$variables)) {
        VARfctText <- paste(VARfctText,"    ",model$variables[[k]]$name," = ",model$variables[[k]]$formula,"\n",sep="")
      }
    }

    if (length(model$reactions)>0) {
      for (k in 1:length(model$reactions)) {
        VARfctText <- paste(VARfctText,"    ",model$reactions[[k]]$name," = ",model$reactions[[k]]$formula,"\n",sep="")
      }
    }

    VARfctText <- paste(VARfctText,"    simresVAR <- rbind(simresVAR,c(\n",sep="")

    if (length(model$variables) > 0) {
      for (k in 1:(length(model$variables))) {
        VARfctText <- paste(VARfctText,"      ",model$variables[[k]]$name,"=unname(",model$variables[[k]]$name,"),\n",sep="")
      }
    }
    if (length(model$reactions) > 0) {
      for (k in 1:(length(model$reactions))) {
        VARfctText <- paste(VARfctText,"      ",model$reactions[[k]]$name,"=unname(",model$reactions[[k]]$name,"),\n",sep="")
      }
    }
    VARfctText <- substr(VARfctText,1,nchar(VARfctText)-2)
    VARfctText <- paste(VARfctText,"))\n",sep="")
    VARfctText <- paste(VARfctText,"  }\n",sep="")
  }

  VARfctText <- paste(VARfctText,"  return(simresVAR)\n",sep="")
  VARfctText <- paste(VARfctText,"}\n",sep="")
  outFunc <- eval(parse(text=VARfctText))
}


###############################################################################
# genROOTsim: Generate event root finding function
###############################################################################
genROOTsim <- function (model) {

  ROOTfctText <- "function (time,states,paramvalues) {\n"
  ROOTfctText <- paste(ROOTfctText,getODEtext(model),sep="")
  ROOTfctText <- paste(ROOTfctText,"  roots <- c(\n",sep="")
  for (k in 1:(getNumberOfEventsAZRmodel(model))) {
    # Check first if trigger function contains either gtAZR, geAZR, ltAZR, leAZR if not throw an error
    trigger <- model$events[[k]]$trigger
    testStart <- c(AZRaux:::strlocateall(trigger,"AZRsim:::le(")$start,
                   AZRaux:::strlocateall(trigger,"AZRsim:::lt(")$start,
                   AZRaux:::strlocateall(trigger,"AZRsim:::ge(")$start,
                   AZRaux:::strlocateall(trigger,"AZRsim:::gt(")$start)
    if (is.null(testStart))
      stop("genROOTsim: Trigger function(s) for event(s) wrongly defined. Please use syntax: gt/ge/lt/le(expr1,expr2)")

    ROOTfctText <- paste(ROOTfctText,"      ",model$events[[k]]$trigger,"-0.5,\n",sep="")
  }
  ROOTfctText <- substr(ROOTfctText,1,nchar(ROOTfctText)-2)
  ROOTfctText <- paste(ROOTfctText,")\n",sep="")
  ROOTfctText <- paste(ROOTfctText,"  return(roots)\n",sep="")
  ROOTfctText <- paste(ROOTfctText,"}\n",sep="")

  outFunc <- eval(parse(text=ROOTfctText))
}


###############################################################################
# genEVASsim: Generate event assignment function
###############################################################################
genEVASsim <- function (model) {

  EVASfctText <- "function (EVENTindex,time,states,paramvalues) {\n"
  EVASfctText <- paste(EVASfctText,getODEtext(model),sep="")

  EVASfctText <- paste(EVASfctText,"  eventAssignment <- c()\n",sep="")
  for (k in 1:getNumberOfEventsAZRmodel(model)) {
    if (getNumberOfEventassignmentsAZRmodel(model,k) == 0)
      stop("genEVASsim: model contains an event without assignments")

    # We generate named vectors, depending on the event index
    EVASfctText <- paste(EVASfctText,"  if (EVENTindex==",k,") {\n",sep="")
    EVASfctText <- paste(EVASfctText,"    eventAssignment <- c(\n",sep="")
    for (k2 in 1:getNumberOfEventassignmentsAZRmodel(model,k)) {
      EVASfctText <- paste(EVASfctText,"      ",
                           model$events[[k]]$assignment[[k2]]$variable," = ",
                           model$events[[k]]$assignment[[k2]]$formula,",\n",sep="")
    }
    EVASfctText <- substr(EVASfctText,1,nchar(EVASfctText)-2)
    EVASfctText <- paste(EVASfctText,")\n",sep="")
    EVASfctText <- paste(EVASfctText,"  }\n",sep="")
  }
  EVASfctText <- paste(EVASfctText,"  return(eventAssignment)\n",sep="")
  EVASfctText <- paste(EVASfctText,"}\n",sep="")
  outFunc <- eval(parse(text=EVASfctText))
}


###############################################################################
# implementALLinputMath: Calls implementInputMath for all inputs in a model
# Returns updated model with all inputs handled.
###############################################################################
implementALLinputMath <- function(model) {
  if (getNumberOfInputsAZRmodel(model)>0) {
    # Implement the math for each input
    for (k in 1:getNumberOfInputsAZRmodel(model)) {
      model <- implementInputMath(model,k)
    }
    # Remove each input from model
    for (k in 1:getNumberOfInputsAZRmodel(model)) {
      # index always 1!
      model <- delInputAZRmodel(model,1)
    }
  }
  return(model)
}

###############################################################################
# implementInputMath: Generates needed model equations and components for a defined INPUT
# This function is called before the creation of the simulation functions
# This function is called from AZRsimulate - but ONLY if a dosing table is provided.
# The dosing table is pre-processed before application!
# How to do this when wanting to estimate some of the parameters still needs to be
# determined - maybe the typical input parameters always can be assumed to be for estimation.
# Returns updated model with input handled.
# NOTE: this function should NOT BE CALLED ALONE but always via implementALLinputMath,
# as it does not remove the input definitions in the model.
###############################################################################
implementInputMath <- function(model,inputindex) {

  # Get input information
  name <- model$inputs[[inputindex]]$name
  factors <- model$inputs[[inputindex]]$factors
  stateindex <- model$inputs[[inputindex]]$stateindex

  # Define parameter names to add for input
  paramDoseAmountName <- paste(name,"dose",sep="")
  paramDoseTimeName <- paste(name,"time",sep="")
  paramDoseDurationName <- paste(name,"duration",sep="")
  paramDoseTlagName <- paste(name,"lagtime",sep="")
  # Add parameters to model with default values
  model <- addParameterAZRmodel(model,paramDoseAmountName,0)
  model <- addParameterAZRmodel(model,paramDoseTimeName,0)
  model <- addParameterAZRmodel(model,paramDoseDurationName,1e-10) # To avoid division by zero
  model <- addParameterAZRmodel(model,paramDoseTlagName,0)

  # Define timing and rate variable
  varDoseName <- tolower(name)
  varDoseFormulaTimingRate <- paste(
      paramDoseAmountName, "/", paramDoseDurationName, "*",
      "piecewise(1,",
        "and(",
          "ge(time,",paramDoseTimeName,"+",paramDoseTlagName,"),",
          "lt(time,",paramDoseTimeName,"+",paramDoseTlagName,"+",paramDoseDurationName,"))",
      ",0)",sep="")
  # Add timing and rate variable to the model
  model <- addVariableAZRmodel(model,varDoseName,formula=varDoseFormulaTimingRate)

  # Add distribution information to the ODEs
  for (k in 1:length(stateindex)) {
    # Add variable to state
    ODE <- model$states[[stateindex[k]]]$ODE
    if (substr(factors[k],1,1) =="+")
      factors[k] <- substr(factors[k],2,nchar(factors[k]))

    if (factors[k]=="1" || factors[k]=="(+1)") {
      ODE <- paste(ODE,"+",varDoseName,sep="")
    } else {
      ODE <- paste(ODE,"+",factors[k],"*",varDoseName,sep="")
    }
    model$states[[stateindex[k]]]$ODE <- ODE
  }

  return(model)
}

