% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagram_helpers.R
\name{prepare_inputs}
\alias{prepare_inputs}
\title{prepare inputs from the system of equations for diagrammer}
\usage{
prepare_inputs(model_states)
}
\arguments{
\item{model_states}{named vector of ode system of equations}
}
\description{
prepare inputs from the system of equations for diagrammer
}
\details{
nodes for diagrammer need to be specified, this will extract the
compartments and the inputs and generate the required syntax
}
