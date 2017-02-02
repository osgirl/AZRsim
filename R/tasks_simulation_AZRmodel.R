###############################################################################
###############################################################################
# Main simulation interface functions for AZRmodels, dosing schemes, etc. Covers
# simulation both in deSolve and CVODES.
###############################################################################
###############################################################################

###############################################################################
# AZRsimulate: Simulate an AZRmodel
###############################################################################
#' Simulate an AZRmodel
#'
#' To be updated! Many things could be added to this simulation.
#' The simulation function expects the model to be augmented with
#' R or C functions for model. If these are not present, the function will add
#' them. If they are present, the function assumes that these are correct and
#' match the model structure.
#'
#' @param model An AZRmodel
#' @param simtime Simulation time vector. If scalar provided then 1001 simulation
#'        steps will be used.
#' @param IC Vector with initial conditions. If not provided the initial
#'        conditions stored in the model will be used as defaults.
#' @param parameters A named parameter vector to be used for simulation.
#'        Named parameters need to be present in the model.
#' @param paramnames An alternative to providing a named vector with parameter
#'        values is to provide a vector with parameter names and a vector
#'        with corresponding parameter values (see next input argument)
#' @param paramvalues A vector with parameter values is 'paramnames' is defined.
#' @param dosingTable A dataframe defining a dosing table. Required columns:
#'        TIME, INPUT, DOSE, DURATION, LAGTIME
#' @param FLAGdosOut If TRUE then dosing table information will be added to the
#'        output variable. If FALSE then it will not be added. This flag only
#'        has effect if a dosingTable is given as input argument.
#' @param method String with the deSolve method for ODE integration (default: lsode
#'        as this is a solver that can reasonably deal with state constraints)
#' @param atol Absolute tolerance for integration
#' @param rtol Relative tolerance for integration
#' @return AZRmodel with attached simulation functions (see above)
#' @examples
#' model <- AZRmodel(system.file("examples/NovakTyson.txt", package="AZRsim"))
#' x <- AZRsimulate(model,400)
#' x <- AZRsimulate(model,400,parameters=c(k1=0.5))
#' x <- AZRsimulate(model,0:400,IC=c(Cyclin=0.5))
#' x <- AZRsimulate(model,0:400,IC=c(Cyclin=0.5),parameters=c(k1=0.9))
#' @export

AZRsimulate <- function (model,
                         simtime = NULL,
                         IC = NULL,
                         parameters = NULL,
                         paramnames = NULL,
                         paramvalues = NULL,
                         dosingTable = NULL,
                         FLAGdosOut = FALSE,
                         method = "lsode",
                         atol = 1e-6,
                         rtol = 1e-6
) {

  if (!is.AZRmodel(model))
    stop("AZRsimulate: provided model argument is not an AZRmodel")

  if (getNumberOfStatesAZRmodel(model) == 0)
    stop("AZRsimulate: provided model has no dynamic states")

  if (hasalgebraicAZRmodel(model))
    stop("AZRsimulate: provided model has algebraic states. This is not supported in AZRsim at the moment")

  if (hasfastreactionsAZRmodel(model))
    stop("AZRsimulate: provided model has fast reactions. This is not supported in AZRsim at the moment")

  if (!is.null(parameters) && (!is.null(paramnames) || !is.null(paramvalues)))
    stop("AZRsimulate: parameters can either be provided by a named vector 'parameters' or by combination of 'paramnames' and 'paramvalues'. But not both at the same time")

  if ((!is.null(paramnames) && is.null(paramvalues)) || (is.null(paramnames) && !is.null(paramvalues)))
    stop("AZRsimulate: if 'paramnames' is defined, also 'paramvalues' needs to be defined - and vice versa!")

  if (length(paramnames) != length(paramvalues))
    stop("AZRsimulate: 'paramnames' and 'paramvalues' need to have same number of elements")

  if (!is.null(dosingTable)) {
    # Check and process dosing table
    dosingTable <- checkProcessDosingTable(dosingTable)
    # Use 1.2 times the max dosing time if a dosing table is defined
    if (is.null(simtime))
      simtime <- max(dosingTable$TIME) * 1.5
  } else {
    # Use 20 time units as default simulation time if no dosing table is defined
    if (is.null(simtime))
        simtime <- 20
  }

  # Handle scalar definition of simtime
  if (length(simtime) == 1)
    simtime <- seq(0,simtime,simtime/1000)

  # Check simtime
  if (length(simtime) != length(unique(simtime)))
    stop("AZRsimulate: simtime vector does not contain unique elements")

  # Check if simulation function are attached to the model
  # If not attached at this point the model was probably loaded with simFlag=FALSE
  # So here we generate the missing information
  if (is.null(attr(model,"ODEsim"))) {
    model <- genSimFunctions(implementALLinputMath(model))
  }

  # Get default parameter values
  paramInfo <- getAllParametersAZRmodel(model)
  parametersDefault <- paramInfo$paramvalues
  names(parametersDefault) <- paramInfo$paramnames

  # Get provided parameters in same format (paramnames, paramvalues)
  if (!is.null(parameters)) {
    paramnames <- names(parameters)
    paramvalues <- parameters
  }

  # Need to check if provided parameter names are all available in the model
  if (!is.null(paramnames)) {
    check <- setdiff(paramnames,names(parametersDefault))
    if (length(check) != 0)
      stop("AZRsimulate: provided parameter names contain names that are not present in the model")
  }

  # Determine simulation parameters
  parametersSim <- parametersDefault
  if (!is.null(paramvalues))
    parametersSim[paramnames] <- paramvalues

  # Determine default ICs and handle potentially non-numerical ICs
  defaultIC <- attr(model,"nnICsim")()

  # Need to check if provided IC names are all available in the model
  if (!is.null(IC)) {
    check <- setdiff(names(IC),names(defaultIC))
    if (length(check) != 0)
      stop("AZRsimulate: provided IC names contain names that are not present in the model")
  }

  # Determine simulation ICs
  ICsim <- defaultIC
  if (!is.null(IC))
    ICsim[names(IC)] <- IC[names(IC)]

  # Here the real simulation starts and different handlings need to be done
  # in the case of availability or non-availability of dosingTable information
  if (is.null(dosingTable)) {
    # Handle case where no dosing table is provided
    simresALL <- simulateAZRmodelDeSolve(model,simtime,ICsim,parametersSim,method,atol,rtol)
  } else {
    simresALL <- simulateAZRmodelDosingTableDeSolve(model,simtime,ICsim,parametersSim,dosingTable,method,atol,rtol)
  }

  # Update the simulation results "time" -> "TIME"
  names(simresALL)[1] <- "TIME"

  # On demand, integrate the dosing table into the simulation results
  # And remove added "inputn" variables
  if (FLAGdosOut && !is.null(dosingTable)) {
    xe <- simresALL
    xe$EVID<-0
    dte <- dosingTable
    dte$EVID<-1
    y <- dplyr::full_join(xe,dte,by=c("TIME","EVID"))
    y <- dplyr::arrange(y,y[,"TIME"])
    y <- dplyr::filter(y,y[,"TIME"]<=max(simtime))
    y$EVID <- NULL
    if (getNumberOfInputsAZRmodel(attr(model,"originalModel")) > 0) {
      for (k in 1:getNumberOfInputsAZRmodel(attr(model,"originalModel"))) {
        y[,paste("input",k,sep="")] <- NULL
      }
    }
    simresALL <- y
  }

  return(simresALL)
}


###############################################################################
# simulateAZRmodelDosingTableDeSolve: Handle dosing table simulation
###############################################################################
simulateAZRmodelDosingTableDeSolve <- function(model,simtime,ICsim,parametersSim,dosingTable,method,atol,rtol) {
  # Cycle through dosing records in dosing table and simulate each dosing period and collect results

  # initialize simresALL
  simresALL <- c()

  # Adjust dosing table to max TIME as in max simtime
  dosingTable <- dplyr::filter(dosingTable,dosingTable[,"TIME"]<=max(simtime))

  # Add information about actual dose administration start
  dosingTable$TIME_DOSE_EFFECT_START <- dosingTable$TIME+dosingTable$LAGTIME

  # Sort by start of dose effect
  dosingTable <- dplyr::arrange(dosingTable,dosingTable[,"TIME_DOSE_EFFECT_START"])

  # Simulations need to be done from current TIME_DOSE_EFFECT_START to next TIME_DOSE_EFFECT_START ...

  # Get unique dosing start time points - can be empty if simulation time too short
  dosingEffectStartTimes <- unique(dosingTable$TIME_DOSE_EFFECT_START)

  # If dosingEffectStartTimes empty then only do first piece in normal way
  if (length(dosingEffectStartTimes)==0) {
    simresALL <- simulateAZRmodelDeSolve(model,simtime,ICsim,parametersSim,method,atol,rtol)
    return(simresALL)
  }

  # Create simulation time vector until and including first dose time
  simtimePreFirstDose <- unique(c(simtime[simtime<dosingEffectStartTimes[1]],dosingEffectStartTimes[1]))

  # Simulate until first dose
  if (length(simtimePreFirstDose)>0) {
    simresPreFirstDose <- simulateAZRmodelDeSolve(model,simtimePreFirstDose,ICsim,parametersSim,method,atol,rtol)
    # Get last state as next initial condition (time point of next dose)
    ICsim <- unlist(simresPreFirstDose[nrow(simresPreFirstDose),2:(getNumberOfStatesAZRmodel(model)+1)])
    # Store simulation results
    # Do not exclude last time point - which is the time of the next dose
    simresALL <- rbind(simresALL,simresPreFirstDose[1:(nrow(simresPreFirstDose)),])
  }

  # Simulate each dose if more than one dose
  if (length(dosingEffectStartTimes) > 1) {
    for (k in 1:(length(dosingEffectStartTimes)-1)) {
      # Create simulation time vector for piece
      # Add as first time point the time for dose and as last the time for next dose
      simtimePiece <- unique(c(dosingEffectStartTimes[k],simtime[simtime>=dosingEffectStartTimes[k] & simtime<=dosingEffectStartTimes[k+1]],dosingEffectStartTimes[k+1]))

      # Get dosing information for the dosing time
      doseInfo <- dplyr::filter(dosingTable,dosingTable[,"TIME_DOSE_EFFECT_START"]==dosingEffectStartTimes[k])

      # Need to generate an updated parameter vector with dosing information
      for (k2 in 1:nrow(doseInfo)) {
        parametersSim[paste("INPUT",doseInfo$INPUT[k2],"dose",sep="")] <- doseInfo$DOSE[k2]
        parametersSim[paste("INPUT",doseInfo$INPUT[k2],"time",sep="")] <- doseInfo$TIME[k2]
        parametersSim[paste("INPUT",doseInfo$INPUT[k2],"duration",sep="")] <- doseInfo$DURATION[k2]
        parametersSim[paste("INPUT",doseInfo$INPUT[k2],"lagtime",sep="")] <- doseInfo$LAGTIME[k2]
      }

      # Simulate piece
      simresPiece <- simulateAZRmodelDeSolve(model,simtimePiece,ICsim,parametersSim,method,atol,rtol)
      # Get last state as next initial condition (time point of next dose)
      ICsim <- unlist(simresPiece[nrow(simresPiece),2:(getNumberOfStatesAZRmodel(model)+1)])
      # Store simulation results
      # Previous piece contained as last entry the dose time and this piece contained as first entry the
      # same dose time. We keep the results from previous piece and remove the first from this piece.
      simresALL <- rbind(simresALL,simresPiece[2:(nrow(simresPiece)),])
    }
  }

  # Simulate last dose until final time
  # Only needed if more than More than 1 time points to simulate remain.

  # Create simulation time vector for time post last dose
  simresPostLastDose <- unique(c(dosingEffectStartTimes[length(dosingEffectStartTimes)],simtime[simtime>=dosingEffectStartTimes[length(dosingEffectStartTimes)]]))

  if (length(simresPostLastDose)>1) {
    # Get dosing information for the dosing time
    doseInfo <- dplyr::filter(dosingTable,dosingTable[,"TIME_DOSE_EFFECT_START"]==dosingEffectStartTimes[length(dosingEffectStartTimes)])

    # Need to generate an updated parameter vector with dosing information
    for (k2 in 1:nrow(doseInfo)) {
      parametersSim[paste("INPUT",doseInfo$INPUT[k2],"dose",sep="")] <- doseInfo$DOSE[k2]
      parametersSim[paste("INPUT",doseInfo$INPUT[k2],"time",sep="")] <- doseInfo$TIME[k2]
      parametersSim[paste("INPUT",doseInfo$INPUT[k2],"duration",sep="")] <- doseInfo$DURATION[k2]
      parametersSim[paste("INPUT",doseInfo$INPUT[k2],"lagtime",sep="")] <- doseInfo$LAGTIME[k2]
    }

    # Simulate piece
    simresPostLastPiece <- simulateAZRmodelDeSolve(model,simresPostLastDose,ICsim,parametersSim,method,atol,rtol)
    # Store simulation results
    # Previous piece contained as last entry the dose time and this piece contained as first entry the
    # same dose time. We keep the results from previous piece and remove the first from this piece.
    simresALL <- rbind(simresALL,simresPostLastPiece[2:(nrow(simresPostLastPiece)),])
  }

  # Keep only simtime elements
  simresALL <- dplyr::filter(simresALL,simresALL[,"time"] %in% simtime)

  return(simresALL)
}


###############################################################################
# simulateAZRmodelDeSolve: Simulate an AZRmodel with deSolve - able to handle
# events ... auxiliary function used with and without dosingTable presence
###############################################################################
simulateAZRmodelDeSolve <- function(model,simtime,ICsim,parametersSim,method,atol,rtol) {
  # Handle simulation for events and no events present in a different manner
  if (getNumberOfEventsAZRmodel(model)==0) {
    # Simulate in deSolve without events
    simresODE <- deSolve::ode(y = ICsim, times = simtime, func = attr(model,"ODEsim"), parms = parametersSim, method=method, atol=atol, rtol=rtol)
    # Determine variable and reaction values
    varresODE <- attr(model,"VARsim")(simresODE,parametersSim)
    # Combine all results to obtain the output
    simresALL <- as.data.frame(cbind(simresODE,varresODE))
  } else {
    # Use an Event-wrapper function to simulate events
    simresALL <- c()
    timeEndPiece <- simtime[1]
    simtimePiece <- simtime
    ICpiece <- ICsim
    parametersSimPiece <- parametersSim
    while(timeEndPiece<max(simtime)) {
      simresODEpiece <- deSolve::ode(y = ICpiece, times = simtimePiece,
                                     func = attr(model,"ODEsim"),
                                     parms = parametersSimPiece,
                                     method=method,
                                     rootfun = attr(model,"ROOTsim"))

      # Determine variable and reaction values for pieces
      varresODEpiece <- attr(model,"VARsim")(simresODEpiece,parametersSimPiece)
      # Combine all results
      simresALLpiece <- as.data.frame(cbind(simresODEpiece,varresODEpiece))
      # Combine result pieces
      simresALL = rbind(simresALL,simresALLpiece)

      # Check if an event has fired
      ixevent <- veclocate(attributes(simresODEpiece)$iroot==1)
      if (length(ixevent!=0)) {
        # Event has fired - get the relevant event assignment function
        # and apply it to states and parameters
        eventAssignmentFct <- attr(model,"EVASsim")
        timeEval <- simresODEpiece[nrow(simresODEpiece),"time"]
        statesEval <- simresODEpiece[nrow(simresODEpiece),2:ncol(simresODEpiece)]
        newValues <- eventAssignmentFct(ixevent,unname(timeEval),unname(statesEval),unname(parametersSimPiece))
        # aux vector
        allElements <- c(statesEval,parametersSimPiece)
        allElements[names(newValues)] <- newValues[names(newValues)]
        # split aux vector
        ICpiece[names(ICpiece)] <- allElements[names(ICpiece)]
        parametersSimPiece[names(parametersSimPiece)] <- allElements[names(parametersSimPiece)]
      }

      # Get new end time from previous step
      timeEndPiece <- simresODEpiece[nrow(simresODEpiece),"time"]
      # Get new simtimePiece vector
      simtimePiece <- c(timeEndPiece,simtimePiece[simtimePiece>timeEndPiece])
    }
    simresALL <- dplyr::filter(simresALL,simresALL[,"time"] %in% simtime)
  }
  return(simresALL)
}
