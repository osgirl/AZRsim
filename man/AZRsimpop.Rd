% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{AZRsimpop}
\alias{AZRsimpop}
\title{Population simulation}
\usage{
AZRsimpop(model, simtime = NULL, icTable = NULL, parameterTable = NULL,
  dosing_table = NULL, FLAGdosOut = FALSE, outputs = NULL,
  FLAGaddParam = FALSE, FLAGaddIC = FALSE, ncores = 1,
  opt_method_stiff = TRUE, opt_abstol = 1e-06, opt_reltol = 1e-06,
  opt_minstep = 0, opt_maxstep = 0, opt_initstep = 0,
  opt_maxnumsteps = 1e+05, opt_maxerrtestfails = 50,
  opt_maxorder_stiff = 5, opt_maxorder_nonstiff = 12,
  opt_maxconvfails = 10, opt_maxnonlineariter = 3, verbose = FALSE)
}
\arguments{
\item{model}{An AZRmodel}

\item{simtime}{Simulation time vector. If scalar provided then 1001 simulation
steps will be used. If not provided (NULL) seq(0,20,1000) will be used if
no dosing_table provided. If not provided and a dosing_table is provided,
seq(0,1.5x the max dosing time,1000) is used. For each individual the
same simulation time vector will be used.}

\item{icTable}{Can be NULL, a vector, a matrix, or a dataframe. In any case the
elements need to be named and each row needs to contain initial condition
information for all states in the model. If NULL, then the ICs stored in
the model will be used. If a vector then for all subjects the same initial
conditions will be used. Otherwise if a matrix or dataframe then each
row corresponds to one subject. If individual initial conditions are provided
then the number of them needs to match potentially individual parameter
and dosing table numbers. Ordering of initial conditions for each subject
needs to be as for parameters and dosing tables. The ID column in the dosing
table is only used for defining individual subjects, not for ordering.
If initial conditions are provided they do override any potential
settings through non-numerical initial conditions in the model.}

\item{parameterTable}{Can be NULL, a vector, a matrix, or a dataframe. In any case the
elements need to be named and need to contain parameters that are present
in the model. Not all parameters need to be defined. Undefined ones are
kept on the value stored in the model. If a vector then for all subjects
the same parameters will be used. Otherwise if a matrix or dataframe then
each row corresponds to one subject. If individual parameters are provided
then the number of them needs to match potentially individual initial
conditions and dosing table numbers. Ordering of parameters for each subject
needs to be as for ICs and dosing tables. The ID column in the dosing
table is only used for defining individual subjects, not for ordering.}

\item{dosing_table}{Can be NULL or a dataframe. If NULL, then no dosing is
simulated. If a datafram is provided it can contain an ID column to
define individual subjects dosings. If no ID column provided, same dosing
used for all subjects. Each dosing is defined as the dosing_table in AZRsimulate.
If individual dosing tables are provided then the number of them needs
to match potentially individual initial conditions and parameter numbers.
Ordering of dosing tables for each subject needs to be as for ICs and
parameters. The ID column in the dosing table is only used for defining
individual subjects, not for ordering.}

\item{FLAGdosOut}{If TRUE then dosing table information will be added to the
output variable. If FALSE then it will not be added. This flag only
has effect if a dosing_table is given as input argument.}

\item{outputs}{A vector with names of outputs to return from simulation.
By default (NULL) all states, variables, reactions, are returned.}

\item{FLAGaddParam}{A flag to indicate if simulation parameters should be added
(TRUE) or not (FALSE) to the output. Only modified / provided parameters
are added.}

\item{FLAGaddIC}{A flag to indicate if simulation ICs should be added
(TRUE) or not (FALSE) to the output.}

\item{ncores}{Number of cores to distribute the individual simulations on}

\item{opt_method_stiff}{Flag (FALSE: non-stiff, TRUE: stiff)}

\item{opt_abstol}{Double value for absolute tolerance}

\item{opt_reltol}{Double value for relative tolerance}

\item{opt_minstep}{Double value for minimal integrator step-size}

\item{opt_maxstep}{Double value for maximal integrator step-size}

\item{opt_initstep}{Double value for initial step-size to be attempted}

\item{opt_maxnumsteps}{Integer value for maximum number of steps between two outputs}

\item{opt_maxerrtestfails}{Integer value for maximum number of error test failures in one step}

\item{opt_maxorder_stiff}{Integer value for maximum order of linear multistep method for STIFF solver (BDF)}

\item{opt_maxorder_nonstiff}{Integer value for maximum order of linear multistep method for NONSTIFF solver (Adams)}

\item{opt_maxconvfails}{Integer value for maximum number of nonlinear solver convergence failures in one step}

\item{opt_maxnonlineariter}{Integer value for maximum number of nonlinear solver iterations permitted per step}

\item{verbose}{Integer flag for outputting additional diagnostic information}
}
\value{
Dataframe with simulation results
}
\description{
Wrapper for the AZRsimulate function, allowing a population simulation in which
individual initial conditions, parameters, and dosing tables can be provided.
Simulation time vector is always the same across simulated subjects. Results
are returned as a dataframe with an ID column, specifying the subjects ID.
}
\examples{
model <- create_model(system.file("examples/NovakTyson.txt", package="AZRsim"))
x <- simulate(model,400)
x <- simulate(model,400,parameters=c(k1=0.5))
}
