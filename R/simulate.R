###############################################################################
###############################################################################
# Main simulation interface functions for AZRmodels, dosing schemes, etc.
# Simulation is done only in CVODES.
###############################################################################
###############################################################################

###############################################################################
# AZRsimulate: Simulate an AZRmodel
###############################################################################
#' Simulate an AZRmodel
#'
#' Simulation function for \code{azrmod} objects which is able to handle dosing events.
#'
#' @param model An object of class \code{azrmod} created using \code{AZRsim::create_model}
#'
#' @param simtime Simulation time vector. If scalar provided then 1001 simulation
#'        steps will be used. If not provided (20) seq(0,20,1000) will be used if
#'        no dosing_table provided. If not provided and a dosing_table is provided,
#'        seq(0,1.5x the max dosing time,1000) is used.
#'
#' @param IC Named vector with numeric initial conditions for ALL states. If not
#'        provided the initial conditions stored in the model will be used as
#'        defaults. Important: If initial conditions are provided, then these
#'        need to be NUMERIC and be provided for ALL states! If the model contains
#'        non-numeric initial conditions, these will be ignored. Basically, user
#'        provided initial conditions overwrite all other settings. If the model
#'        itself does contain non-numeric initial conditions, then it might be
#'        more useful to change these via parameter settings, using the "parameters"
#'        input argument.
#'
#' @param parameters A named parameter vector to be used for simulation.
#'        Named parameters need to be present in the model. Initial condition
#'        definitions in AZRmodels can be non-numeric mathematical expressions
#'        and depend on states, parameters, variables, reactions, and functions.
#'        A states IC can depend on a previously defined state in the model.
#'        Otherwise an error message will appear during the import of a model from
#'        text.
#'
#' @param dosing_table A dataframe defining a dosing table with the following columns:
#'        \code{TIME}, \code{INPUT}, \code{DOSE}, \code{DURATION}, and \code{LAGTIME}.
#' @param FLAGdosOut If TRUE then dosing table information will be added to the
#'        output variable. If FALSE then it will not be added. This flag only
#'        has effect if a dosing_table is given as input argument.
#'
#' @param outputs A vector with names of outputs to return from simulation.
#'        By default (NULL) all states, variables, reactions, are returned.
#'
#' @param opt_method_stiff      Flag (FALSE: non-stiff, TRUE: stiff)
#' @param opt_abstol            Double value for absolute tolerance
#' @param opt_reltol            Double value for relative tolerance
#' @param opt_minstep           Double value for minimal integrator step-size
#' @param opt_maxstep           Double value for maximal integrator step-size
#' @param opt_initstep          Double value for initial step-size to be attempted
#' @param opt_maxnumsteps       Integer value for maximum number of steps between two outputs
#' @param opt_maxerrtestfails   Integer value for maximum number of error test failures in one step
#' @param opt_maxorder_stiff    Integer value for maximum order of linear multistep method for STIFF solver (BDF)
#' @param opt_maxorder_nonstiff Integer value for maximum order of linear multistep method for NONSTIFF solver (Adams)
#' @param opt_maxconvfails      Integer value for maximum number of nonlinear solver convergence failures in one step
#' @param opt_maxnonlineariter  Integer value for maximum number of nonlinear solver iterations permitted per step
#' @param verbose               Integer flag for outputting additional diagnostic information
#'
#' @return A \code{data.frame} object of class \code{azrsim} which contains simulations.
#' @examples
#' # simple harmonic oscillator simulation
#' sho_model <- create_model(system.file("examples/sho.txt", package="AZRsim"))
#' sho_sim <- simulate(sho_model, seq(1, 100, by = 0.1))
#' sho_sim <- simulate(sho_model, 100, parameters = c("theta" = 0.5))
#' sho_sim <- simulate(sho_model, 100, IC = c("y1" = 1, "y2" = 2))
#' plot(sho_sim)
#'
#' # simple one compartment dosing
#' one_cpt <- create_model(system.file("examples/one_cpt_dt.txt", package="AZRsim"))
#' dt <- data.frame("TIME" = seq(1,9, by = 1),
#'                  "DOSE" = 40,
#'                  "DURATION" = 0,
#'                  "INPUT" = 1,
#'                  "LAGTIME" = 0,
#'                  stringsAsFactors = FALSE)
#' one_cpt_sim <- simulate(one_cpt, seq(0, 10, by=0.01), dosing_table = dt, output = c("y"))
#' plot(one_cpt_sim, lwd = 2, plot_names = "blood")
#' @export

simulate.azrmod <- function (model,
                             # Simulation time vector
                             simtime               = NULL,
                             # Initial conditions
                             IC                    = NULL,
                             # Parameter information
                             parameters            = NULL,
                             # Dosing / Event information
                             dosing_table           = NULL,
                             FLAGdosOut            = FALSE,
                             # Output definitions
                             outputs               = NULL,
                             # Define integrator etc. options
                             opt_method_stiff      = TRUE,
                             opt_abstol            = 1.0e-6,
                             opt_reltol            = 1.0e-6,
                             opt_minstep           = 0.0,
                             opt_maxstep           = 0.0,
                             opt_initstep          = 0.0,
                             opt_maxnumsteps       = 100000,
                             opt_maxerrtestfails   = 50,
                             opt_maxorder_stiff    = 5,
                             opt_maxorder_nonstiff = 12,
                             opt_maxconvfails      = 10,
                             opt_maxnonlineariter  = 3,
                             verbose               = FALSE
) {

  ##############################################################################
  # Basic AZRmodel checks
  ##############################################################################

  if (!is_azrmod(model))
    stop("AZRsimulate: provided model argument is not an AZRmodel")

  if (len_states(model) == 0)
    stop("AZRsimulate: provided model has no dynamic states")

  if (has_algebraic(model))
    stop("AZRsimulate: provided model has algebraic states. This is not supported in AZRsim at the moment")

  if (has_fast_reactions(model))
    stop("AZRsimulate: provided model has fast reactions. This is not supported in AZRsim at the moment")

  ##############################################################################
  # Check if simulation function is attached to the model
  # If not attached at this point the model was probably loaded with simFlag=FALSE
  ##############################################################################

  if (is.null(attr(model,"modelCfunAddress"))) model <- genSimFunctions(model)

  ##############################################################################
  # Check if model DLL is still loaded  otherwise load it again
  ##############################################################################
  if (!attr(model,"modelDLLname") %in% names(getLoadedDLLs()))
    stop("AZRsimulate: Something happened - model DLL got unloaded. Please reload the model")

  ##############################################################################
  # Check outputs input argument
  ##############################################################################

  if (!is.null(outputs)) {
    test <- setdiff(outputs,c(get_all_states(model)$statenames,get_all_variables(model)$varnames,get_all_reactions(model)$reacnames))
    if (length(test) != 0)
      stop("AZRsimulate: At least one element defined in 'outputs' is not present in the model as state, variable, or reaction")
  }

  ##############################################################################
  # Check and handle simtime settings
  ##############################################################################

  if (is.null(simtime)) {
    # Use 20 time units as default simulation time if no dosing table is defined
    simtime <- 20
    # and 1.5x the max dosing time if a dosing_table is defined
    if (!is.null(dosing_table)) {
      max_dose_time <- max(dosing_table$TIME)
      if (is.null(max_dose_time)) stop("AZRsimulate: dosing_table provided without TIME column")
      simtime <- max_dose_time * 1.5
      if (max_dose_time==0) simtime <- 20
    }
  }
  # If simtime is a scalar then use 1001 equidistant time points for simulation
  if (length(simtime) == 1) simtime <- seq(0,simtime,simtime/1000)
  # Final check of simtime vector
  if (length(simtime) != length(unique(simtime))) stop("AZRsimulate: simtime vector does not contain unique elements")

  ##############################################################################
  # Check and process dosing_table information - if undefined it remains undefined
  ##############################################################################

  dosing_table <- check_dosing_table(dosing_table)

  ##############################################################################
  # Handle simulation parameter values
  ##############################################################################

  # Get default parameter values stored in the AZRmodel
  parameters_default <- get_all_parameters(model)$paramvalues

  # Need to check if provided parameter names are all available in the model
  if (!is.null(parameters)) {
    check <- setdiff(names(parameters),names(parameters_default))
    if (length(check) != 0) stop("AZRsimulate: provided parameter names contain names that are not present in the model")
  }

  # Determine simulation parameters based on default parameters and user requested parameters
  parameters_sim <- parameters_default
  if (!is.null(parameters)) parameters_sim[names(parameters)] <- parameters

  ##############################################################################
  # Handle simulation initial conditions
  ##############################################################################

  # Get default initial conditions, taking into account potentially non-numeric
  # initial conditions that might depend on parameters and varibles.
  defaultIC <- calcNNic(model,parameters_sim)

  # If initial conditions are provided by the user then ensure that they are provided for ALL states
  # And bring them into the right order of states (IC needs to be a named vector)
  # User provided initial conditions will overwrite any other initial condition settings.
  if (!is.null(IC)) {
    # Check that ICs provided for all states in the model
    if (length(IC) != length(defaultIC))
      stop("AZRsimulate: when providing initial conditions, values for all states need to be provided as a named vector")
    if (length(setdiff(names(defaultIC),names(IC))) != 0)
      stop("AZRsimulate: provided IC names contain names that are not present in the model")
    if (!isnumericVector(IC))
      stop("AZRsimulate: User-provided initial condition vector needs to be numeric.")
    # Assign new values and take care of order via the names
    ICsim <- defaultIC
    ICsim[names(IC)] <- IC
  } else {
    # If not user provided then use the default initial conditions
    # parameter settings will have been taken into account in case of non-numerical initial conditions
    ICsim <- defaultIC
  }

  ##############################################################################
  # Determine needed input arguments for the CVODES C interface
  ##############################################################################

  # Get component number information
  NRSTATES      <- len_states(model)
  NRPARAMETERS  <- len_parameters(model)
  NRVARIABLES   <- len_variables(model)
  NRREACTIONS   <- len_reactions(model)
  NREVENTS      <- len_events(model)
  model_elements_nr <- c(NRSTATES,NRPARAMETERS,NRVARIABLES,NRREACTIONS,NREVENTS)

  # Get C code models address
  model_func_ptr <- attr(model,"modelCfunAddress")

  ##############################################################################
  # Call simulation functions - without or with dosing_table
  ##############################################################################

  if (is.null(dosing_table)) {
    # Handle case where no dosing table is provided
    # Call cvodes integrator interface
    simres_all <- .Call("cvodesAZRinterface",             # Name of C-code CVODES interface function
                       PACKAGE="AZRsim",                    # Name of the DLL file in which the interface function is located
                       model_func_ptr,                      # Pointer to model function
                       as.integer(model_elements_nr),       # Integer vector with numbers of model elements
                       as.double(simtime),                  # Double vector with time points for simulation
                       as.double(ICsim),                    # Double vector with initial conditions
                       as.double(parameters_sim),            # Double vector with parameter values
                       as.integer(opt_method_stiff),        # Integer flag (0: non-stiff, 1:stiff)
                       as.double(opt_abstol),               # Double value for absolute tolerance
                       as.double(opt_reltol),               # Double value for relative tolerance
                       as.double(opt_minstep),              # Double value for minimal integrator step-size
                       as.double(opt_maxstep),              # Double value for maximal integrator step-size
                       as.double(opt_initstep),             # Double value for initial step-size to be attempted
                       as.integer(opt_maxnumsteps),         # Integer value for maximum number of steps between two outputs
                       as.integer(opt_maxerrtestfails),     # Integer value for maximum number of error test failures in one step
                       as.integer(opt_maxorder_stiff),      # Integer value for maximum order of linear multistep method for STIFF solver (BDF)
                       as.integer(opt_maxorder_nonstiff),   # Integer value for maximum order of linear multistep method for NONSTIFF solver (Adams)
                       as.integer(opt_maxconvfails),        # Integer value for maximum number of nonlinear solver convergence failures in one step
                       as.integer(opt_maxnonlineariter),    # Integer value for maximum number of nonlinear solver iterations permitted per step
                       as.integer(FALSE),                   # Integer value defining what to do: 0=do integration, 1=return RHS of ODE for given time[0], states, parameters
                       as.integer(verbose)                  # Integer flag for outputting additional diagnostic information
    )
  } else {
    simres_all <- simulateAZRmodelDosingTable(model_func_ptr,
                                             model_elements_nr,
                                             simtime,
                                             ICsim,
                                             parameters_sim,
                                             dosing_table,
                                             opt_method_stiff,
                                             opt_abstol,
                                             opt_reltol,
                                             opt_minstep,
                                             opt_maxstep,
                                             opt_initstep,
                                             opt_maxnumsteps,
                                             opt_maxerrtestfails,
                                             opt_maxorder_stiff,
                                             opt_maxorder_nonstiff,
                                             opt_maxconvfails,
                                             opt_maxnonlineariter,
                                             verbose)
  }

  ##############################################################################
  # Post-process simulation output to named dataframe
  ##############################################################################

  # Convert simulation output to dataframe
  simres_all <- as.data.frame(simres_all)

  # Update names
  names(simres_all) <- c("TIME", get_all_states(model)$statenames, get_all_variables(model)$varnames,get_all_reactions(model)$reacnames)

  # Keep only simtime elements
  simres_all <- dplyr::filter(simres_all,TIME %in% simtime)

  #############################################################################
  # Handle outputs if defined
  #############################################################################

  if (!is.null(outputs)) {
    # Keep only TIME and elements in outputs
    simres_all <- simres_all[,c("TIME",outputs)]
  }

  #############################################################################
  # Post-process simulation output by adding dose information if desired by user
  #############################################################################
  # On demand, integrate the dosing table into the simulation results
  # And remove added "inputn" variables
  if (FLAGdosOut && !is.null(dosing_table)) {
    xe          <- simres_all
    xe$EVID     <- 0
    dte         <- dosing_table
    dte$EVID    <- 1
    y           <- dplyr::full_join(xe,dte,by=c("TIME","EVID"))
    y           <- dplyr::arrange(y,TIME)
    y           <- dplyr::filter(y,TIME<=max(simtime))
    y$EVID      <- NULL
    if (len_inputs(attr(model,"originalModel")) > 0) {
      for (k in 1:len_inputs(attr(model,"originalModel"))) {
        y[,paste("input",k,sep="")] <- NULL
      }
    }
    simres_all    <- y
  }
  class(simres_all) <- c("azrsim", "data.frame")
  return(simres_all)
}


###############################################################################
# simulateAZRmodelDosingTable: Handle dosing table simulation
# This function should be included in the C interface in a future optimization
# of the code. It pbly would speed up simulation of dosing scenarios a little.
###############################################################################
simulateAZRmodelDosingTable <- function(model_func_ptr,
                                        model_elements_nr,
                                        simtime,
                                        ICsim,
                                        parameters_sim,
                                        dosing_table,
                                        opt_method_stiff,
                                        opt_abstol,
                                        opt_reltol,
                                        opt_minstep,
                                        opt_maxstep,
                                        opt_initstep,
                                        opt_maxnumsteps,
                                        opt_maxerrtestfails,
                                        opt_maxorder_stiff,
                                        opt_maxorder_nonstiff,
                                        opt_maxconvfails,
                                        opt_maxnonlineariter,
                                        verbose) {
  # Cycle through dosing records in dosing table and simulate each dosing period and collect results

  # Get number of states
  NRSTATES <- model_elements_nr[1]

  # initialize simres_all
  simres_all <- c()

  # Adjust dosing table to max TIME as in max simtime
  dosing_table <- dplyr::filter(dosing_table,dosing_table[,"TIME"]<=max(simtime))

  # Add information about actual dose administration start
  dosing_table$TIME_DOSE_EFFECT_START <- dosing_table$TIME+dosing_table$LAGTIME

  # Sort by start of dose effect
  dosing_table <- dplyr::arrange(dosing_table,dosing_table[,"TIME_DOSE_EFFECT_START"])

  # Simulations need to be done from current TIME_DOSE_EFFECT_START to next TIME_DOSE_EFFECT_START ...

  # Get unique dosing start time points - can be empty if simulation time too short
  dosing_effect_start_times <- unique(dosing_table$TIME_DOSE_EFFECT_START)

  # If dosing_effect_start_times empty then only do first piece in normal way
  if (length(dosing_effect_start_times)==0) {
    simres_all <- .Call("cvodesAZRinterface",                # Name of C-code CVODES interface function
                       PACKAGE="AZRsim",                    # Name of the DLL file in which the interface function is located
                       model_func_ptr,                      # Pointer to model function
                       as.integer(model_elements_nr),       # Integer vector with numbers of model elements
                       as.double(simtime),                  # Double vector with time points for simulation
                       as.double(ICsim),                    # Double vector with initial conditions
                       as.double(parameters_sim),            # Double vector with parameter values
                       as.integer(opt_method_stiff),        # Integer flag (0: non-stiff, 1:stiff)
                       as.double(opt_abstol),               # Double value for absolute tolerance
                       as.double(opt_reltol),               # Double value for relative tolerance
                       as.double(opt_minstep),              # Double value for minimal integrator step-size
                       as.double(opt_maxstep),              # Double value for maximal integrator step-size
                       as.double(opt_initstep),             # Double value for initial step-size to be attempted
                       as.integer(opt_maxnumsteps),         # Integer value for maximum number of steps between two outputs
                       as.integer(opt_maxerrtestfails),     # Integer value for maximum number of error test failures in one step
                       as.integer(opt_maxorder_stiff),      # Integer value for maximum order of linear multistep method for STIFF solver (BDF)
                       as.integer(opt_maxorder_nonstiff),   # Integer value for maximum order of linear multistep method for NONSTIFF solver (Adams)
                       as.integer(opt_maxconvfails),        # Integer value for maximum number of nonlinear solver convergence failures in one step
                       as.integer(opt_maxnonlineariter),    # Integer value for maximum number of nonlinear solver iterations permitted per step
                       as.integer(FALSE),                   # Integer value defining what to do: 0=do integration, 1=return RHS of ODE for given time[0], states, parameters
                       as.integer(verbose)                  # Integer flag for outputting additional diagnostic information
    )
    return(simres_all)
  }

  # Create simulation time vector until and including first dose time
  simtimePreFirstDose <- unique(c(simtime[simtime<dosing_effect_start_times[1]],
                                  dosing_effect_start_times[1]))

  # Simulate until first dose
  if (length(simtimePreFirstDose)>1) {
    simresPreFirstDose <- .Call("cvodesAZRinterface",       # Name of C-code CVODES interface function
                                PACKAGE="AZRsim",                    # Name of the DLL file in which the interface function is located
                                model_func_ptr,                      # Pointer to model function
                                as.integer(model_elements_nr),       # Integer vector with numbers of model elements
                                as.double(simtimePreFirstDose),      # Double vector with time points for simulation
                                as.double(ICsim),                    # Double vector with initial conditions
                                as.double(parameters_sim),            # Double vector with parameter values
                                as.integer(opt_method_stiff),        # Integer flag (0: non-stiff, 1:stiff)
                                as.double(opt_abstol),               # Double value for absolute tolerance
                                as.double(opt_reltol),               # Double value for relative tolerance
                                as.double(opt_minstep),              # Double value for minimal integrator step-size
                                as.double(opt_maxstep),              # Double value for maximal integrator step-size
                                as.double(opt_initstep),             # Double value for initial step-size to be attempted
                                as.integer(opt_maxnumsteps),         # Integer value for maximum number of steps between two outputs
                                as.integer(opt_maxerrtestfails),     # Integer value for maximum number of error test failures in one step
                                as.integer(opt_maxorder_stiff),      # Integer value for maximum order of linear multistep method for STIFF solver (BDF)
                                as.integer(opt_maxorder_nonstiff),   # Integer value for maximum order of linear multistep method for NONSTIFF solver (Adams)
                                as.integer(opt_maxconvfails),        # Integer value for maximum number of nonlinear solver convergence failures in one step
                                as.integer(opt_maxnonlineariter),    # Integer value for maximum number of nonlinear solver iterations permitted per step
                                as.integer(FALSE),                   # Integer value defining what to do: 0=do integration, 1=return RHS of ODE for given time[0], states, parameters
                                as.integer(verbose)                  # Integer flag for outputting additional diagnostic information
    )
    # Get last state as next initial condition (time point of next dose)
    ICsim <- unlist(simresPreFirstDose[nrow(simresPreFirstDose),2:(NRSTATES+1)])
    # Store simulation results
    # Do not exclude last time point - which is the time of the next dose
    simres_all <- rbind(simres_all,simresPreFirstDose[1:(nrow(simresPreFirstDose)),])
    addFirst = FALSE
  } else {
    addFirst = TRUE
  }

  # Simulate each dose if more than one dose
  if (length(dosing_effect_start_times) > 1) {
    for (k in 1:(length(dosing_effect_start_times)-1)) {
      # Create simulation time vector for piece
      # Add as first time point the time for dose and as last the time for next dose
      simtime_piece <- unique(c(dosing_effect_start_times[k],
                               simtime[simtime>=dosing_effect_start_times[k] & simtime<=dosing_effect_start_times[k+1]],
                               dosing_effect_start_times[k+1]))

      # Get dosing information for the dosing time
      doseInfo <- dplyr::filter(dosing_table,dosing_table[,"TIME_DOSE_EFFECT_START"]==dosing_effect_start_times[k])

      # Need to generate an updated parameter vector with dosing information
      for (k2 in 1:nrow(doseInfo)) {
        parameters_sim[paste("INPUT",doseInfo$INPUT[k2],"dose",sep="")] <- doseInfo$DOSE[k2]
        parameters_sim[paste("INPUT",doseInfo$INPUT[k2],"time",sep="")] <- doseInfo$TIME[k2]
        parameters_sim[paste("INPUT",doseInfo$INPUT[k2],"duration",sep="")] <- doseInfo$DURATION[k2]
        parameters_sim[paste("INPUT",doseInfo$INPUT[k2],"lagtime",sep="")] <- doseInfo$LAGTIME[k2]
      }

      # Simulate piece
      simres_piece <- .Call("cvodesAZRinterface",                       # Name of C-code CVODES interface function
                           PACKAGE="AZRsim",                    # Name of the DLL file in which the interface function is located
                           model_func_ptr,                      # Pointer to model function
                           as.integer(model_elements_nr),       # Integer vector with numbers of model elements
                           as.double(simtime_piece),             # Double vector with time points for simulation
                           as.double(ICsim),                    # Double vector with initial conditions
                           as.double(parameters_sim),            # Double vector with parameter values
                           as.integer(opt_method_stiff),        # Integer flag (0: non-stiff, 1:stiff)
                           as.double(opt_abstol),               # Double value for absolute tolerance
                           as.double(opt_reltol),               # Double value for relative tolerance
                           as.double(opt_minstep),              # Double value for minimal integrator step-size
                           as.double(opt_maxstep),              # Double value for maximal integrator step-size
                           as.double(opt_initstep),             # Double value for initial step-size to be attempted
                           as.integer(opt_maxnumsteps),         # Integer value for maximum number of steps between two outputs
                           as.integer(opt_maxerrtestfails),     # Integer value for maximum number of error test failures in one step
                           as.integer(opt_maxorder_stiff),      # Integer value for maximum order of linear multistep method for STIFF solver (BDF)
                           as.integer(opt_maxorder_nonstiff),   # Integer value for maximum order of linear multistep method for NONSTIFF solver (Adams)
                           as.integer(opt_maxconvfails),        # Integer value for maximum number of nonlinear solver convergence failures in one step
                           as.integer(opt_maxnonlineariter),    # Integer value for maximum number of nonlinear solver iterations permitted per step
                           as.integer(FALSE),                   # Integer value defining what to do: 0=do integration, 1=return RHS of ODE for given time[0], states, parameters
                           as.integer(verbose)                  # Integer flag for outputting additional diagnostic information
      )

      # Get last state as next initial condition (time point of next dose)
      ICsim <- unlist(simres_piece[nrow(simres_piece),2:(NRSTATES+1)])
      # Store simulation results
      # Previous piece contained as last entry the dose time and this piece contained as first entry the
      # same dose time. We keep the results from previous piece and remove the first from this piece.
      if (addFirst) {
        simres_all <- rbind(simres_all,simres_piece[1:(nrow(simres_piece)),])
      } else {
        simres_all <- rbind(simres_all,simres_piece[2:(nrow(simres_piece)),])
      }

      addFirst = FALSE
    }
  }

  # Simulate last dose until final time
  # Only needed if more than More than 1 time points to simulate remain.

  # Create simulation time vector for time post last dose
  simtime_post_last_dose <- unique(c(dosing_effect_start_times[length(dosing_effect_start_times)],
                                  simtime[simtime>=dosing_effect_start_times[length(dosing_effect_start_times)]]))

  if (length(simtime_post_last_dose)>1) {
    # Get dosing information for the dosing time
    doseInfo <- dplyr::filter(dosing_table,dosing_table[,"TIME_DOSE_EFFECT_START"]==dosing_effect_start_times[length(dosing_effect_start_times)])

    # Need to generate an updated parameter vector with dosing information
    for (k2 in 1:nrow(doseInfo)) {
      parameters_sim[paste("INPUT",doseInfo$INPUT[k2],"dose",sep="")] <- doseInfo$DOSE[k2]
      parameters_sim[paste("INPUT",doseInfo$INPUT[k2],"time",sep="")] <- doseInfo$TIME[k2]
      parameters_sim[paste("INPUT",doseInfo$INPUT[k2],"duration",sep="")] <- doseInfo$DURATION[k2]
      parameters_sim[paste("INPUT",doseInfo$INPUT[k2],"lagtime",sep="")] <- doseInfo$LAGTIME[k2]
    }

    # Simulate piece
    simres_post_last_piece <- .Call("cvodesAZRinterface",        # Name of C-code CVODES interface function
                                 PACKAGE="AZRsim",                    # Name of the DLL file in which the interface function is located
                                 model_func_ptr,                      # Pointer to model function
                                 as.integer(model_elements_nr),       # Integer vector with numbers of model elements
                                 as.double(simtime_post_last_dose),      # Double vector with time points for simulation
                                 as.double(ICsim),                    # Double vector with initial conditions
                                 as.double(parameters_sim),            # Double vector with parameter values
                                 as.integer(opt_method_stiff),        # Integer flag (0: non-stiff, 1:stiff)
                                 as.double(opt_abstol),               # Double value for absolute tolerance
                                 as.double(opt_reltol),               # Double value for relative tolerance
                                 as.double(opt_minstep),              # Double value for minimal integrator step-size
                                 as.double(opt_maxstep),              # Double value for maximal integrator step-size
                                 as.double(opt_initstep),             # Double value for initial step-size to be attempted
                                 as.integer(opt_maxnumsteps),         # Integer value for maximum number of steps between two outputs
                                 as.integer(opt_maxerrtestfails),     # Integer value for maximum number of error test failures in one step
                                 as.integer(opt_maxorder_stiff),      # Integer value for maximum order of linear multistep method for STIFF solver (BDF)
                                 as.integer(opt_maxorder_nonstiff),   # Integer value for maximum order of linear multistep method for NONSTIFF solver (Adams)
                                 as.integer(opt_maxconvfails),        # Integer value for maximum number of nonlinear solver convergence failures in one step
                                 as.integer(opt_maxnonlineariter),    # Integer value for maximum number of nonlinear solver iterations permitted per step
                                 as.integer(FALSE),                   # Integer value defining what to do: 0=do integration, 1=return RHS of ODE for given time[0], states, parameters
                                 as.integer(verbose)                  # Integer flag for outputting additional diagnostic information
    )

    # Store simulation results
    # Previous piece contained as last entry the dose time and this piece contained as first entry the
    # same dose time. We keep the results from previous piece and remove the first from this piece.

    if (addFirst) {
      simres_all <- rbind(simres_all,simres_post_last_piece[1:(nrow(simres_post_last_piece)),])
    } else {
      simres_all <- rbind(simres_all,simres_post_last_piece[2:(nrow(simres_post_last_piece)),])
    }

  }

  return(simres_all)
}


###############################################################################
# AZRxdotcalc: Return evaluated RHS of AZRmodel ODEs
###############################################################################
#' Return evaluated RHS of AZRmodel ODEs
#'
#' Calculates and returns the RHS of ODEs in an AZRmodel for given states,
#' parameters, and time (needs to be scalar). Dosing tables not allowed.
#' Useful to assess issues in a larger AZRmodel.
#'
#' @param model An AZRmodel
#'
#' @param time  Scalar time at which to evaluate the RHS
#'
#' @param states Names vector with state values at which to evaluate the ODE RHS
#'        Values for all states in the model need to be provided. If not provided
#'        then the initial conditions in the model will be used. If non-numeric then
#'        potential changes through parameter settings are taken into account.
#'
#' @param parameters A named parameter vector to be used for evaluation
#'        Named parameters need to be present in the model.
#'
#' @return Vector with evaluated RHS of ODEs
#' @examples
#' model <- AZRmodel(system.file("examples/NovakTyson.txt", package="AZRsim"))
#' x <- AZRxdotcalc(model)
#' x <- AZRxdotcalc(model,400)
#' x <- AZRxdotcalc(model,400,parameters=c(k1=0.5))
#' @export

AZRxdotcalc <- function (model,
                         time                  = 0,
                         states                = NULL,
                         parameters            = NULL) {

  ##############################################################################
  # Basic AZRmodel checks
  ##############################################################################

  if (!is_azrmod(model))
    stop("AZRxdotcalc: provided model argument is not an AZRmodel")

  if (!len_states(model))
    stop("AZRxdotcalc: provided model has no dynamic states")

  if (has_algebraic(model))
    stop("AZRxdotcalc: provided model has algebraic states. This is not supported in AZRsim at the moment")

  if (has_fast_reactions(model))
    stop("AZRxdotcalc: provided model has fast reactions. This is not supported in AZRsim at the moment")

  ##############################################################################
  # Check if simulation function is attached to the model
  # If not attached at this point the model was probably loaded with simFlag=FALSE
  ##############################################################################

  if (is.null(attr(model,"modelCfunAddress"))) model <- genSimFunctions(model)

  ##############################################################################
  # Handle scalar definition of time
  ##############################################################################

  if (length(time) != 1)
    stop("AZRxdotcalc: time needs to be a scalar")

  ##############################################################################
  # Handle simulation parameter values
  ##############################################################################

  # Get default parameter values stored in the AZRmodel
  parameters_default        <- get_all_parameters(model)$paramvalues

  # Need to check if provided parameter names are all available in the model
  if (!is.null(parameters)) {
    check <- setdiff(names(parameters),names(parameters_default))
    if (length(check) != 0) stop("AZRxdotcalc: provided parameter names contain names that are not present in the model")
  }

  # Determine simulation parameters based on default parameters and user requested parameters
  parameters_sim <- parameters_default
  if (!is.null(parameters)) parameters_sim[names(parameters)] <- parameters

  ##############################################################################
  # Handle simulation initial conditions
  ##############################################################################

  # Get default states for Xdot calculation, taking into account potentially non-numeric
  # initial conditions that might depend on parameters and varibles.
  defaultStatesXdotCalc <- calcNNic(model,parameters_sim)

  # If state values are provided by the user then ensure that they are provided for ALL states
  # And bring them into the right order of states (states needs to be a named vector)
  # User provided state values will overwrite any other state settings.
  if (!is.null(states)) {
    # Check that states provided for all states in the model
    if (length(states) != length(defaultStatesXdotCalc))
      stop("AZRxdotcalc: when providing state values, values for all states need to be provided as a named vector")
    if (length(setdiff(names(defaultStatesXdotCalc),names(states))) != 0)
      stop("AZRxdotcalc: provided state names contain names that are not present in the model")
    if (!isnumericVector(states))
      stop("AZRxdotcalc: User-provided states vector needs to be numeric.")
    # Assign new values and take care of order via the names
    statesXdotCalc <- defaultStatesXdotCalc
    statesXdotCalc[names(states)] <- states
  } else {
    # If not user provided then use the default initial conditions
    # parameter settings will have been taken into account in case of non-numerical initial conditions
    statesXdotCalc <- defaultStatesXdotCalc
  }

  ##############################################################################
  # Determine needed input arguments for the CVODES C interface
  ##############################################################################

  # Get component number information
  NRSTATES      <- len_states(model)
  NRPARAMETERS  <- len_parameters(model)
  NRVARIABLES   <- len_variables(model)
  NRREACTIONS   <- len_reactions(model)
  NREVENTS      <- len_events(model)
  model_elements_nr <- c(NRSTATES,NRPARAMETERS,NRVARIABLES,NRREACTIONS,NREVENTS)

  # Get C code models address
  model_func_ptr <- attr(model,"modelCfunAddress")

  ##############################################################################
  # Call cvodes integrator interface
  ##############################################################################

  simres_all <- .Call("cvodesAZRinterface",                # Name of C-code CVODES interface function
                     PACKAGE="AZRsim",                    # Name of the DLL file in which the interface function is located
                     model_func_ptr,                      # Pointer to model function
                     as.integer(model_elements_nr),       # Integer vector with numbers of model elements
                     as.double(time),                     # Double vector with time points for simulation
                     as.double(statesXdotCalc),           # Double vector with initial conditions
                     as.double(parameters_sim),            # Double vector with parameter values
                     as.integer(TRUE),                    # Integer flag (0: non-stiff, 1:stiff)
                     as.double(1.0e-6),                   # Double value for absolute tolerance
                     as.double(1.0e-6),                   # Double value for relative tolerance
                     as.double(0.0),                      # Double value for minimal integrator step-size
                     as.double(0.0),                      # Double value for maximal integrator step-size
                     as.double(0.0),                      # Double value for initial step-size to be attempted
                     as.integer(100000),                  # Integer value for maximum number of steps between two outputs
                     as.integer(50),                      # Integer value for maximum number of error test failures in one step
                     as.integer(5),                       # Integer value for maximum order of linear multistep method for STIFF solver (BDF)
                     as.integer(12),                      # Integer value for maximum order of linear multistep method for NONSTIFF solver (Adams)
                     as.integer(10),                      # Integer value for maximum number of nonlinear solver convergence failures in one step
                     as.integer(3),                       # Integer value for maximum number of nonlinear solver iterations permitted per step
                     as.integer(TRUE),                    # Integer value defining what to do: 0=do integration, 1=return RHS of ODE for given time[0], states, parameters
                     as.integer(FALSE)                    # Integer flag for outputting additional diagnostic information
  )

  # Update names
  names(simres_all) <- c(get_all_states(model)$statenames)

  # Return results
  return(simres_all)
}


##############################################################################
# Evaluate potentially non-numerical initial conditions subject to potential
# parameter changes
##############################################################################

calcNNic <- function(model,parameters_sim) {

  # Check if model contains non-numerical initial conditions ... if not just return the
  # numerical ones
  if (has_only_numeric_ic(model)) return(get_all_states(model)$stateICs)

  # Model contains non-numerical initial conditions => evaluate them - taking into account
  # potential changes in the parameters
  calcIC <- attr(model,"nnICsim")(parameters_sim)

  # Return result
  return(calcIC)
}


###############################################################################
# AZRsimpop: Population simulation
###############################################################################
#' Population simulation
#'
#' Wrapper for the AZRsimulate function, allowing a population simulation in which
#' individual initial conditions, parameters, and dosing tables can be provided.
#' Simulation time vector is always the same across simulated subjects. Results
#' are returned as a dataframe with an ID column, specifying the subjects ID.
#'
#' @param model An AZRmodel
#'
#' @param simtime Simulation time vector. If scalar provided then 1001 simulation
#'        steps will be used. If not provided (NULL) seq(0,20,1000) will be used if
#'        no dosing_table provided. If not provided and a dosing_table is provided,
#'        seq(0,1.5x the max dosing time,1000) is used. For each individual the
#'        same simulation time vector will be used.
#'
#' @param icTable Can be NULL, a vector, a matrix, or a dataframe. In any case the
#'        elements need to be named and each row needs to contain initial condition
#'        information for all states in the model. If NULL, then the ICs stored in
#'        the model will be used. If a vector then for all subjects the same initial
#'        conditions will be used. Otherwise if a matrix or dataframe then each
#'        row corresponds to one subject. If individual initial conditions are provided
#'        then the number of them needs to match potentially individual parameter
#'        and dosing table numbers. Ordering of initial conditions for each subject
#'        needs to be as for parameters and dosing tables. The ID column in the dosing
#'        table is only used for defining individual subjects, not for ordering.
#'        If initial conditions are provided they do override any potential
#'        settings through non-numerical initial conditions in the model.
#'
#' @param parameterTable Can be NULL, a vector, a matrix, or a dataframe. In any case the
#'        elements need to be named and need to contain parameters that are present
#'        in the model. Not all parameters need to be defined. Undefined ones are
#'        kept on the value stored in the model. If a vector then for all subjects
#'        the same parameters will be used. Otherwise if a matrix or dataframe then
#'        each row corresponds to one subject. If individual parameters are provided
#'        then the number of them needs to match potentially individual initial
#'        conditions and dosing table numbers. Ordering of parameters for each subject
#'        needs to be as for ICs and dosing tables. The ID column in the dosing
#'        table is only used for defining individual subjects, not for ordering.
#'
#' @param dosing_table Can be NULL or a dataframe. If NULL, then no dosing is
#'        simulated. If a datafram is provided it can contain an ID column to
#'        define individual subjects dosings. If no ID column provided, same dosing
#'        used for all subjects. Each dosing is defined as the dosing_table in AZRsimulate.
#'        If individual dosing tables are provided then the number of them needs
#'        to match potentially individual initial conditions and parameter numbers.
#'        Ordering of dosing tables for each subject needs to be as for ICs and
#'        parameters. The ID column in the dosing table is only used for defining
#'        individual subjects, not for ordering.
#'
#' @param FLAGdosOut If TRUE then dosing table information will be added to the
#'        output variable. If FALSE then it will not be added. This flag only
#'        has effect if a dosing_table is given as input argument.
#'
#' @param outputs A vector with names of outputs to return from simulation.
#'        By default (NULL) all states, variables, reactions, are returned.
#'
#' @param FLAGaddParam A flag to indicate if simulation parameters should be added
#'        (TRUE) or not (FALSE) to the output. Only modified / provided parameters
#'        are added.
#'
#' @param FLAGaddIC A flag to indicate if simulation ICs should be added
#'        (TRUE) or not (FALSE) to the output.
#'
#' @param ncores Number of cores to distribute the individual simulations on
#'
#' @param opt_method_stiff      Flag (FALSE: non-stiff, TRUE: stiff)
#' @param opt_abstol            Double value for absolute tolerance
#' @param opt_reltol            Double value for relative tolerance
#' @param opt_minstep           Double value for minimal integrator step-size
#' @param opt_maxstep           Double value for maximal integrator step-size
#' @param opt_initstep          Double value for initial step-size to be attempted
#' @param opt_maxnumsteps       Integer value for maximum number of steps between two outputs
#' @param opt_maxerrtestfails   Integer value for maximum number of error test failures in one step
#' @param opt_maxorder_stiff    Integer value for maximum order of linear multistep method for STIFF solver (BDF)
#' @param opt_maxorder_nonstiff Integer value for maximum order of linear multistep method for NONSTIFF solver (Adams)
#' @param opt_maxconvfails      Integer value for maximum number of nonlinear solver convergence failures in one step
#' @param opt_maxnonlineariter  Integer value for maximum number of nonlinear solver iterations permitted per step
#' @param verbose               Integer flag for outputting additional diagnostic information
#'
#' @return Dataframe with simulation results
#' @examples
#' model <- AZRmodel(system.file("examples/NovakTyson.txt", package="AZRsim"))
#' x <- AZRsimulate(model,400)
#' x <- AZRsimulate(model,400,parameters=c(k1=0.5))
#' @export

AZRsimpop <- function (model,
                       # Simulation time vector
                       simtime               = NULL,
                       # Initial conditions
                       icTable               = NULL,
                       # Parameter information
                       parameterTable        = NULL,
                       # Dosing / Event information
                       dosing_table           = NULL,
                       FLAGdosOut            = FALSE,
                       # Outputs
                       outputs               = NULL,
                       FLAGaddParam          = FALSE,
                       FLAGaddIC             = FALSE,
                       # Parallelization settings
                       ncores                = 1,
                       # Define integrator etc. options
                       opt_method_stiff      = TRUE,
                       opt_abstol            = 1.0e-6,
                       opt_reltol            = 1.0e-6,
                       opt_minstep           = 0.0,
                       opt_maxstep           = 0.0,
                       opt_initstep          = 0.0,
                       opt_maxnumsteps       = 100000,
                       opt_maxerrtestfails   = 50,
                       opt_maxorder_stiff    = 5,
                       opt_maxorder_nonstiff = 12,
                       opt_maxconvfails      = 10,
                       opt_maxnonlineariter  = 3,
                       verbose               = FALSE
) {

  ##############################################################################
  # Basic AZRmodel checks
  ##############################################################################

  if (!is_azrmod(model))
    stop("AZRsimpop: provided model argument is not an AZRmodel")

  if (len_states(model) == 0)
    stop("AZRsimpop: provided model has no dynamic states")

  if (has_algebraic(model))
    stop("AZRsimpop: provided model has algebraic states. This is not supported in AZRsim at the moment")

  if (has_fast_reactions(model))
    stop("AZRsimpop: provided model has fast reactions. This is not supported in AZRsim at the moment")

  ##############################################################################
  # Check if simulation function is attached to the model
  # If not attached at this point the model was probably loaded with simFlag=FALSE
  ##############################################################################

  if (is.null(attr(model,"modelCfunAddress"))) model <- genSimFunctions(model)

  ##############################################################################
  # Check and handle simtime settings
  ##############################################################################

  if (is.null(simtime)) {
    # Use 20 time units as default simulation time if no dosing table is defined
    simtime <- 20
    # and 1.5x the max dosing time if a dosing_table is defined
    if (!is.null(dosing_table)) {
      max_dose_time <- max(dosing_table$TIME)
      if (is.null(max_dose_time)) stop("AZRsimpop: dosing_table provided without TIME column")
      simtime <- max_dose_time * 1.5
    }
  }
  # If simtime is a scalar then use 1001 equidistant time points for simulation
  if (length(simtime) == 1) simtime <- seq(0,simtime,simtime/1000)
  # Final check of simtime vector
  if (length(simtime) != length(unique(simtime))) stop("AZRsimulate: simtime vector does not contain unique elements")

  ##############################################################################
  # Handle icTable input argument and bring to dataframe
  ##############################################################################

  # If ICs are provided then ensure they are brought to a dataframe
  if (!is.null(icTable)) {
    if (is.vector(icTable)) {
      namesCol <- names(icTable)
      icTable <- as.data.frame(matrix(icTable,nrow=1))
      colnames(icTable) <- namesCol
    } else {
      if (is.matrix(icTable)) {
        icTable <- as.data.frame(icTable)
      } else {
        if (!is.data.frame(icTable))
          stop("AZRsimpop: type of provided icTable input argument is not correct")
      }
    }
    # Get number ICs
    NR_ICS <- nrow(icTable)
  } else {
    # No icTable provided
    NR_ICS <- 0
  }

  # Check if ID column present ... we do not allow it to not allow people to think that
  # a match based on IDs is done between ICs, parameters, and dosings
  if ("ID" %in% colnames(icTable))
    stop("AZRsimpop: ID column present in icTable. This is not allowed - to avoid thinking ICs, parameters, and dosings are matched by ID!")

  ##############################################################################
  # Handle parameterTable input argument and bring to dataframe
  ##############################################################################

  # If parameters are provided then ensure they are brought to a dataframe
  if (!is.null(parameterTable)) {
    if (is.vector(parameterTable)) {
      namesCol <- names(parameterTable)
      parameterTable <- as.data.frame(matrix(parameterTable,nrow=1))
      colnames(parameterTable) <- namesCol
    } else {
      if (is.matrix(parameterTable)) {
        parameterTable <- as.data.frame(parameterTable)
      } else {
        if (!is.data.frame(parameterTable))
          stop("AZRsimpop: type of provided parameterTable input argument is not correct")
      }
    }
    # Get number parameters
    NR_PARAMETERS <- nrow(parameterTable)
  } else {
    # No parameterTable provided
    NR_PARAMETERS <- 0
  }

  # Check if ID column present ... we do not allow it to not allow people to think that
  # a match based on IDs is done between ICs, parameters, and dosings
  if ("ID" %in% colnames(parameterTable))
    stop("AZRsimpop: ID column present in parameterTable. This is not allowed - to avoid thinking ICs, parameters, and dosings are matched by ID!")

  ##############################################################################
  # Handle/check dosing_table
  # We do not need to handle dosing table completely ... this is done for each
  # subject in the AZRsimulate function ... only sanity checks
  ##############################################################################

  if (!is.null(dosing_table)) {
    # dosing_table provided
    # If given it needs to be a dataframe
    if (!is.data.frame(dosing_table))
      stop("AZRsimpop: provided dosing_table is not a dataframe")

    # Check if ID column present - it is OK not to have an ID column and in this case
    # the same dosing will be given to all subjects
    if ("ID" %in% colnames(dosing_table)) {
      NR_DOSINGTABLE <- length(unique(dosing_table$ID))
    } else {
      NR_DOSINGTABLE <- 1
    }
  } else {
    # No dosing_table given
    NR_DOSINGTABLE <- 0
  }

  ##############################################################################
  # Check if numbers make sense
  ##############################################################################

  # Check first NR_ICS and NR_PARAMETERS and get NR_SUBJECTS info
  if (NR_ICS <= 1) {
    NR_SUBJECTS <- max(NR_PARAMETERS,1)
  } else {
    if (NR_PARAMETERS <= 1) {
      NR_SUBJECTS <- max(NR_ICS,1)
    } else {
      if (NR_ICS != NR_PARAMETERS) stop("AZRsimpop: Number of rows in icTable needs to match number of rows in parameterTable")
      NR_SUBJECTS <- NR_PARAMETERS
    }
  }

  # Now check NR_SUBJECTS against NR_DOSINGTABLE
  if (NR_SUBJECTS == 1) {
    NR_SUBJECTS <- max(NR_DOSINGTABLE,1)
  } else {
    if (NR_DOSINGTABLE > 1) {
      if (NR_DOSINGTABLE != NR_SUBJECTS) stop("AZRsimpop: wrong number of subject level dosing information in dosing_table")
    }
  }

  ##############################################################################
  # Simulation loop
  ##############################################################################

  ##############################################################################
  # Simulation loop - Handle parallelization if desired
  ##############################################################################
  k = NA # avoid CRAN NOTE

  if (ncores==1) {

    "%do%" <- foreach::"%do%"

    # Sequential evaluation
    simresAll <- foreach::foreach (k=1:NR_SUBJECTS, .combine="rbind") %do% {
      indivSimulation(NR_ICS,
                      icTable,
                      NR_PARAMETERS,
                      parameterTable,
                      NR_DOSINGTABLE,
                      dosing_table,
                      k,
                      model,
                      simtime,
                      outputs,
                      FLAGdosOut,
                      opt_method_stiff,
                      opt_abstol,
                      opt_reltol,
                      opt_minstep,
                      opt_maxstep,
                      opt_initstep,
                      opt_maxnumsteps,
                      opt_maxerrtestfails,
                      opt_maxorder_stiff,
                      opt_maxorder_nonstiff,
                      opt_maxconvfails,
                      opt_maxnonlineariter,
                      verbose,
                      FLAGaddParam,
                      FLAGaddIC)
    }
  } else {
    # Parallel execution

    # Set DLL to not loaded
    DLLloaded <- FALSE

    # Start cluster
    maxCores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(ncores,maxCores))
    doParallel::registerDoParallel(cl)

    # Export complete environment to cluster workers
    parallel::clusterExport(cl,envir=environment(),varlist=ls())

    "%dopar%" <- foreach::"%dopar%"

    # Parallel foreach loop
    simresAll <- foreach::foreach (k=1:NR_SUBJECTS, .packages="AZRsim", .combine="rbind", .inorder=TRUE) %dopar% {

      if (!DLLloaded) {
        DLLfile <- attr(model,"modelDLLfile")
        if (.Platform$OS.type=="unix") {
          DLLfile <- paste0(DLLfile,".so")
        }
        dyn.load(DLLfile)
        model_func_ptr <- getNativeSymbolInfo("model",PACKAGE=attr(model,"modelDLLname"))$address
        attr(model,"modelCfunAddress") <- model_func_ptr
        DLLloaded <- TRUE
      }

      indivSimulation(NR_ICS,
                      icTable,
                      NR_PARAMETERS,
                      parameterTable,
                      NR_DOSINGTABLE,
                      dosing_table,
                      k,
                      model,
                      simtime,
                      outputs,
                      FLAGdosOut,
                      opt_method_stiff,
                      opt_abstol,
                      opt_reltol,
                      opt_minstep,
                      opt_maxstep,
                      opt_initstep,
                      opt_maxnumsteps,
                      opt_maxerrtestfails,
                      opt_maxorder_stiff,
                      opt_maxorder_nonstiff,
                      opt_maxconvfails,
                      opt_maxnonlineariter,
                      verbose,
                      FLAGaddParam,
                      FLAGaddIC)
    }
    # Stop cluster
    suppressWarnings(doParallel::stopImplicitCluster())
  }

  # Update rownames sequentially
  rownames(simresAll) <- 1:nrow(simresAll)

  ##############################################################################
  # Return
  ##############################################################################

  return(simresAll)
}


##############################################################################
# Content of individual simulation loop ... own function to allow to switch
# between foreach and for depending on number of cores
##############################################################################
indivSimulation <- function(NR_ICS,
                            icTable,
                            NR_PARAMETERS,
                            parameterTable,
                            NR_DOSINGTABLE,
                            dosing_table,
                            k,
                            model,
                            simtime,
                            outputs,
                            FLAGdosOut,
                            opt_method_stiff,
                            opt_abstol,
                            opt_reltol,
                            opt_minstep,
                            opt_maxstep,
                            opt_initstep,
                            opt_maxnumsteps,
                            opt_maxerrtestfails,
                            opt_maxorder_stiff,
                            opt_maxorder_nonstiff,
                            opt_maxconvfails,
                            opt_maxnonlineariter,
                            verbose,
                            FLAGaddParam,
                            FLAGaddIC) {

  # Determine individual initial conditions
  if (NR_ICS == 0) {
    # None provided
    icIndiv <- NULL
  } else {
    if (NR_ICS == 1) {
      # Single one provided (use for all)
      icIndiv <- as.numeric(icTable[1,])
    } else {
      # Multiple provided (one per subject)
      icIndiv <- as.numeric(icTable[k,])
    }
    # Add column names again
    names(icIndiv) <- colnames(icTable)
    # If provided then potentially available non-numeric initial conditions in the
    # model are overridden!
  }

  # Determine individual parameters
  if (NR_PARAMETERS == 0) {
    parametersIndiv <- NULL
  } else {
    if (NR_PARAMETERS == 1) {
      # Single one provided (use for all)
      parametersIndiv <- as.numeric(parameterTable[1,])
    } else {
      # Multiple provided (one per subject)
      parametersIndiv <- as.numeric(parameterTable[k,])
    }
    # Add column names again
    names(parametersIndiv) <- colnames(parameterTable)
  }

  # Determine individual dosing table
  if (NR_DOSINGTABLE == 0) {
    dosing_tableIndiv <- NULL
  } else {
    if (NR_DOSINGTABLE == 1) {
      dosing_tableIndiv <- dosing_table
    } else {
      # In contrast to MATLAB unique() in R does not sort ... which is nice!
      allID <- unique(dosing_table$ID)
      dosing_tableIndiv <- dplyr::filter(dosing_table,ID==allID[k])
    }
    # Remove ID if present
    dosing_tableIndiv$ID <- NULL
  }

  # Simulate individual
  simresIndiv <- simulate(model                 = model,
                             simtime               = simtime,
                             IC                    = icIndiv,
                             parameters            = parametersIndiv,
                             dosing_table           = dosing_tableIndiv,
                             FLAGdosOut            = FLAGdosOut,
                             outputs               = outputs,
                             opt_method_stiff      = opt_method_stiff,
                             opt_abstol            = opt_abstol,
                             opt_reltol            = opt_reltol,
                             opt_minstep           = opt_minstep,
                             opt_maxstep           = opt_maxstep,
                             opt_initstep          = opt_initstep,
                             opt_maxnumsteps       = opt_maxnumsteps,
                             opt_maxerrtestfails   = opt_maxerrtestfails,
                             opt_maxorder_stiff    = opt_maxorder_stiff,
                             opt_maxorder_nonstiff = opt_maxorder_nonstiff,
                             opt_maxconvfails      = opt_maxconvfails,
                             opt_maxnonlineariter  = opt_maxnonlineariter,
                             verbose               = verbose)

  # Construct individual result
  # Always add ID
  simresIndiv <- cbind(ID=k, simresIndiv)
  # Add initial conditions if desired
  if (FLAGaddIC && !is.null(icTable)) {
    if (NR_ICS==1) {
      simresIndiv <- cbind(simresIndiv, icTable[rep(1,nrow(simresIndiv)),])
    } else {
      simresIndiv <- cbind(simresIndiv, icTable[rep(k,nrow(simresIndiv)),])
    }
  }
  # Add parameters if desired
  if (FLAGaddParam && !is.null(parameterTable)) {
    if (NR_PARAMETERS == 1) {
      simresIndiv <- cbind(simresIndiv, parameterTable[rep(1,nrow(simresIndiv)),])
    } else {
      simresIndiv <- cbind(simresIndiv, parameterTable[rep(k,nrow(simresIndiv)),])
    }
  }

  # Return results
  return(simresIndiv)
}
