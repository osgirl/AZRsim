###############################################################################
###############################################################################
# Generation of needed simulation functions for AZRmodels.
###############################################################################
###############################################################################


###############################################################################
# genSimFunctions: Generate simulation functions for the C simulator interface
#                  and attach them to the model
###############################################################################
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
  # Generate and add non-numerical IC handling function
  ##############################################################################
  # For this potentially available interpolation function need to get some
  # syntax change so that non-numerical initial conditions can be evaluated in R
  # even in the case of presence of interpolations.
  #   Replace syntax for interpolation functions from
  #   interpx([comma sep vector],[comma sep vector],element) to
  #   interpx(c(comma sep vector),c(comma sep vector),element)
  #   More in general all "vector definitions" with "[...]" are replaced by "c(...)"
  modelInterpSyntaxR <- handleVectorSyntaxInterpR(model)
  attr(model,"nnICsim") <- nnICsim(modelInterpSyntaxR)

  ##############################################################################
  # Generate C code model, compile it, load it, attach address of model function
  # to attributes
  ##############################################################################
  model <- compileAZRmodelAZR(model)

  ##############################################################################
  # Return the updated model
  ##############################################################################
  return(model)
}

###############################################################################
# compileAZRmodelAZR: Converts an AZRmodel to C-code and compiles it
###############################################################################
# Converts an AZRmodel to C-code and compiles it
#
# Converts an AZRmodel to C code, compiles it to DLL and loads this DLL.
# Attaches information about DLL path and model address to the AZRmodel object.
# If previous DLL present, it will attempt to unload the previous DLL.
#
# Function normally not needed - except if model object was
#
# @param model An AZRmodel
# @return AZRmodel with attached address and DLL information

compileAZRmodelAZR <- function(model) {

  if (!is.AZRmodel(model))
    stop("compileAZRmodelAZR: input argument is not an AZRmodel")

  # Check if model was already compiled and in this case unload DLL and
  # remove the DLL file
  if (!is.null(attr(model,"modelDLLfile"))) {
    try(dyn.unload(attr(model,"modelDLLfile")), silent=TRUE)
    try(unlink(paste(attr(model,"modelDLLfile"),".*",sep="")),silent=TRUE)
  }

  # Get temporary file name for C code model and DLL
  modelCfilepath <- tempfile()
  modelCpath     <- fileparts(modelCfilepath)$pathname
  modelCfilename <- fileparts(modelCfilepath)$filename

  # Change to temporary folder
  oldpath        <- getwd()
  setwd(modelCpath)

  # Convert to C
  exportCcodeAZRmodel(model,modelCfilename)

  # Create Makevars file in same folder as the model is located (work folder)
  includesPaths <- .libPaths()
  includesLocDef <- paste('PKG_CPPFLAGS =',sep="")
  for (k in seq_along(includesPaths))
    includesLocDef <- paste(includesLocDef,' -I"',includesPaths[k],'/AZRsim/solver/include/"',sep="")
  includesLocDef <- paste(includesLocDef, "\nPKG_LIBS=  -lm")

  readr::write_file(text=includesLocDef,filename="Makevars")

  # Compile model to DLL - stdout to xxx
  system(paste("R CMD SHLIB ",modelCfilename,".c",sep=""))

  # Clean files
  unlink(paste(modelCfilename,".c",sep=""))
  unlink(paste(modelCfilename,".o",sep=""))
  unlink("Makevars")

  # Load model DLL
  tryCatch({
    if (.Platform$OS.type=="unix") {
      dyn.load(paste0(modelCfilename,".so"))
    } else {
      dyn.load(modelCfilename)
    }
  }, error = function(err) {
    setwd(oldpath)
    stop("compileAZRmodelAZR: Compilation of model led to an error. Please check the above output and correct the model.")
  })

  # Return to old working directory
  setwd(oldpath)

  # Get pointer to rhs function of C-code ODE model
  model_func_ptr <- getNativeSymbolInfo("model",PACKAGE=modelCfilename)$address

  # Add DLL information to the model
  attr(model,"modelDLLfile")     <- modelCfilepath
  attr(model,"modelDLLname")     <- modelCfilename
  attr(model,"modelCfunAddress") <- model_func_ptr

  # Return model
  return(model)
}

###############################################################################
# handleVectorSyntaxInterpR: Handles vector syntax
###############################################################################
handleVectorSyntaxInterpR <- function(model) {

  # Handle ODEs
  for (k in 1:getNumberOfStatesAZRmodel(model)) {
    formula <- model$states[[k]]$ODE
    # Check if square brackets present in formula - if so then lets
    # make R vectors out of it
    if (!is.null(strlocateall(formula,"[")$start)) {
      # Square brackets present
      formula <- strremWhite(formula)
      formula <- strrepM(formula,"[","c(")
      formula <- strrepM(formula,"]",")")
      model$states[[k]]$ODE <- formula
    }
  }

  # Handle Variables
  if (getNumberOfVariablesAZRmodel(model) > 0) {
    for (k in 1:getNumberOfVariablesAZRmodel(model)) {
      formula <- model$variables[[k]]$formula
      # Check if square brackets present in formula - if so then lets
      # make R vectors out of it
      if (!is.null(strlocateall(formula,"[")$start)) {
        # Square brackets present
        formula <- strremWhite(formula)
        formula <- strrepM(formula,"[","c(")
        formula <- strrepM(formula,"]",")")
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
      if (!is.null(strlocateall(formula,"[")$start)) {
        # Square brackets present
        formula <- strremWhite(formula)
        formula <- strrepM(formula,"[","c(")
        formula <- strrepM(formula,"]",")")
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

###############################################################################
# nnICsim: Handle non-numeric initial conditions
# This function generates a function for the evaluation of non-numerical initial
# conditions. This generated function is faster to evaluate and the good thing is
# that here during the generation of this function it is already checked if the
# non-numerical initial conditions can be evaluated ... meaning during import of
# a model rather than during simulation of an imported model. Also the possible
# dependencies of nn ICs are much wider than in the IQM Tools.
###############################################################################
nnICsim <- function (model) {

  ###############################################
  # STEP 1: Expand all IC RHSs to only states and parameters
  ###############################################
  # GET ALL FORMULAS AND ODES
  stateInfo <- getAllStatesAZRmodel(model)
  paramInfo <- getAllParametersAZRmodel(model)
  varInfo   <- getAllVariablesAZRmodel(model)
  reacInfo  <- getAllReactionsAZRmodel(model)
  # EXPAND VARIABLES
  if (length(varInfo$varnames) > 0) {
    for (k in 1:length(varInfo$varnames)) {
      if (k+1 <= length(varInfo$varnames)) {
        for (k2 in (k+1):length(varInfo$varnames)) {
          varInfo$varformulas[k2] <- gsub(pattern = paste("\\<",varInfo$varnames[k],"\\>",sep=""),replacement = paste("(",varInfo$varformulas[k],")",sep=""),x = varInfo$varformulas[k2])
        }
      }
    }
  }
  # EXPAND REACTIONS
  if (length(reacInfo$reacnames) > 0) {
    for (k in 1:length(reacInfo$reacnames)) {
      if (k+1 <= length(reacInfo$reacnames)) {
        for (k2 in (k+1):length(reacInfo$reacnames)) {
          reacInfo$reacformulas[k2] <- gsub(pattern = paste("\\<",reacInfo$reacnames[k],"\\>",sep=""),replacement = paste("(",reacInfo$reacformulas[k],")",sep=""),x = reacInfo$reacformulas[k2])
        }
      }
    }
  }
  # INSERT REACTIONS INTO ODES (might contain VARIABLES)
  if (length(stateInfo$statenames) > 0) {
    for (k in 1:length(stateInfo$statenames)) {
      if (length(reacInfo$reacnames) > 0) {
        for (k2 in 1:length(reacInfo$reacnames)) {
          updatedIC <- gsub(pattern = paste("\\<",reacInfo$reacnames[k2],"\\>",sep=""),replacement = paste("(",reacInfo$reacformulas[k2],")",sep=""),x = stateInfo$stateICs[k])
          stateInfo$stateICs[k] <- updatedIC
        }
      }
    }
  }
  # INSERT VARIABLES INTO ODES
  if (length(stateInfo$statenames) > 0) {
    for (k in 1:length(stateInfo$statenames)) {
      if (length(varInfo$varnames) > 0) {
        for (k2 in 1:length(varInfo$varnames)) {
          updatedIC <- gsub(pattern = paste("\\<",varInfo$varnames[k2],"\\>",sep=""),replacement = paste("(",varInfo$varformulas[k2],")",sep=""),x = stateInfo$stateICs[k])
          stateInfo$stateICs[k] <- updatedIC
        }
      }
    }
  }

  ###############################################
  # STEP 2: Test run the equations to see if they can be evaluated ... if not then
  #         throw an error and tell user that non-numeric ICs badly formed
  ###############################################
  # Define functions
  for (k in seq_along(model$functions)) {
    text <- paste(model$functions[[k]]$name," <- function(",model$functions[[k]]$arguments,") { ",model$functions[[k]]$formula," }")
    eval(parse(text=text))
  }
  # Initialize states with NA
  for (k in seq_along(stateInfo$statenames)) {
    text <- paste(stateInfo$statenames[k],"= NA")
    eval(parse(text=text))
  }
  # Initialize parameters with default values
  for (k in seq_along(paramInfo$paramnames)) {
    text <- paste(paramInfo$paramnames[k],"=",paramInfo$paramvalues[k])
    eval(parse(text=text))
  }
  # Evaluate IC RHS expressions in order of appearance in the model
  for (k in seq_along(stateInfo$statenames)) {
    text <- paste(stateInfo$statenames[k],"=",stateInfo$stateICs[k])
    eval(parse(text=text))
  }
  # Collect resulting numerical ICs
  test <- c()
  for (k in seq_along(stateInfo$statenames)) {
    text <- paste("test[k] <- ",stateInfo$statenames[k])
    eval(parse(text=text))
  }
  # Check if any of the elements in "test" is NA ... then NN ICs are non-evaluable
  if (sum(as.double(is.na(test))) > 0)
    stop("nnICsim: Non-numerical initial conditions are wrongly defined and non-evaluable.")

  ###############################################
  # STEP 3: Generate NN ICs function
  ###############################################

  nnICfctText <- "function (parameters) {\n"

  # Define functions
  for (k in seq_along(model$functions))
    nnICfctText <- paste(nnICfctText,"  ",model$functions[[k]]$name," <- function(",model$functions[[k]]$arguments,") { ",model$functions[[k]]$formula," }\n",sep="")

  # Initialize states with NA
  for (k in seq_along(stateInfo$statenames))
    nnICfctText <- paste(nnICfctText,"  ",stateInfo$statenames[k]," <- NA\n",sep="")

  # Initialize parameters with provided values
  for (k in seq_along(paramInfo$paramnames))
    nnICfctText <- paste(nnICfctText,"  ",paramInfo$paramnames[k]," <- parameters[",k,"]\n",sep="")

  # Evaluate IC RHS expressions
  for (k in seq_along(stateInfo$statenames))
    nnICfctText <- paste(nnICfctText,"  ",stateInfo$statenames[k]," <- ",stateInfo$stateICs[k],"\n",sep="")

  # Collect results
  nnICfctText <- paste(nnICfctText,"  out <- c()\n",sep="")
  for (k in seq_along(stateInfo$statenames))
    nnICfctText <- paste(nnICfctText,'  out["',stateInfo$statenames[k],'"] <- ',stateInfo$statenames[k],"\n",sep="")

  # Check if any of the elements in "test" is NA ... then NN ICs are non-evaluable
  # This should not happen except if parameters are wrongly defined ...
  nnICfctText <- paste(nnICfctText,"  if (sum(as.double(is.na(out))) > 0) stop('Non-numerical initial conditions are wrongly defined and non-evaluable.')\n",sep="")

  # Return
  nnICfctText <- paste(nnICfctText,"  return(out)\n",sep="")
  nnICfctText <- paste(nnICfctText,"}\n",sep="")
  outFunc <- eval(parse(text=nnICfctText))

  ###############################################
  # STEP 4: Done! Return the generated function
  ###############################################
  return(outFunc)
}
