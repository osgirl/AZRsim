###############################################################################
###############################################################################
# This file contains the export functions for AZRmodels to TXT files
# These are not exported and called from export_azrmod
###############################################################################
###############################################################################

###############################################################################
# exportTxtAZRmodel: exports AZRmodel as .txt file
###############################################################################
# Export of AZRmodel to .txt file
#
# Export of an AZRmodel to a .txt file. ODE representation is used. These
# .txt files allow a humanly readable format of the model files. An alternative
# is the export to a .txtbc file, which uses a chemical reaction type of
# syntax, if possible. The default filename is constructed from the models name.
#
# @param model An AZRmodel to be exported.
# @param filename Full path with filename to export the model to (.txt will be added)
# @return None

exportTxtAZRmodel <- function (model, filename=NULL) {

  if (!is_azrmod(model))
    stop("exportTxtAZRmodel: input argument is not an AZRmodel")

  if (is.null(filename))
    filename <- gsub("\\W","",model$name)

  filename <- paste(strrepM(filename,".txt",""), ".txt", sep="")

  # Initialize the FILETEXT
  FILETEXT <- ""

  FILETEXT <- paste(FILETEXT,"********** MODEL NAME\n\n",model$name,"\n\n",sep="")

  FILETEXT <- paste(FILETEXT,"********** MODEL NOTES\n\n",model$notes,"\n\n",sep="")

  ###############################################
  # States
  ###############################################

  FILETEXT <- paste(FILETEXT,"********** MODEL STATES\n\n",sep="")

  if (getNumberOfStatesAZRmodel(model) > 0) {
    for (k in 1:getNumberOfStatesAZRmodel(model)) {
      ODEtext <- paste("d/dt(",model$states[[k]]$name,") = ",model$states[[k]]$ODE,sep="")
      type <- model$states[[k]]$type
      compartment <- model$states[[k]]$compartment
      unittype <- model$states[[k]]$unittype
      informationText <- ""

      if (!is.null(type) || !is.null(compartment) || !is.null(unittype)) {
        if (type=="isSpecie" && !is.null(compartment) && (unittype=="amount" || unittype=="concentration"))
          informationText <- paste(" {",type,":",compartment,":",unittype,"}",sep="")
        if (type=="isParameter" && is.null(compartment) && is.null(unittype))
          informationText <- paste(" {",type,"}",sep="")
        if (type=="isCompartment" && is.null(unittype))
          informationText <- paste(" {",type,":",compartment,"}",sep="")
        if (informationText=="") {
          stop(paste("exportTxtAZRmodel: Type information for state ",model$states[[k]]$name," seems to be wrong."),sep="")
        }
      }

      if (!is.null(model$states[[k]]$lowConstraint) && !is.null(model$states[[k]]$highConstraint)) {
        constraintsText = paste(" {constraints:[",model$states[[k]]$lowConstraint,",",model$states[[k]]$highConstraint,"]}",sep="")
      } else {
        constraintsText = ""
      }

      ODEtext <- strtrimM(paste(ODEtext,informationText,constraintsText,sep=""))
      if (!is.null(model$states[[k]]$notes))
        ODEtext <- strtrimM(paste(ODEtext,"%",model$states[[k]]$notes,sep=" "))
      FILETEXT <- paste(FILETEXT,ODEtext,"\n",sep="")
    }
    FILETEXT <- paste(FILETEXT," \n",sep="")
  }

  ###############################################
  # Algebraic Rules
  ###############################################

  for (k in seq_along(model$algebraic)) {
    if (!is.null(model$algebraic[[k]]$name)) {
      ALGtext <- paste("0 = ",model$algebraic[[k]]$formula," : ",model$algebraic[[k]]$name,sep="")
    } else {
      ALGtext <- paste("0 = ",model$algebraic[[k]]$formula,sep="")
    }
    type <- model$algebraic[[k]]$type
    compartment <- model$algebraic[[k]]$compartment
    unittype <- model$algebraic[[k]]$unittype
    informationText <- ""
    if (!is.null(type) || !is.null(compartment) || !is.null(unittype)) {
      if (type=="isSpecie" && !is.null(compartment) && (unittype=="amount" || unittype=="concentration"))
        informationText <- paste(" {",type,":",compartment,":",unittype,"}",sep="")
      if (type=="isParameter" && is.null(compartment) && is.null(unittype))
        informationText <- paste(" {",type,"}",sep="")
      if (type=="isCompartment" && is.null(unittype))
        informationText <- paste(" {",type,":",compartment,"}",sep="")
      if (informationText=="") {
        stop(paste("exportTxtAZRmodel: Type information for algebraic state ",
                   model$algebraic[[k]]$name, " seems to be wrong."), sep="")
      }
    }
    ALGtext <- strtrimM(paste(ALGtext,informationText,sep=""))
    if (!is.null(model$algebraic[[k]]$notes))
      ALGtext <- strtrimM(paste(ALGtext,"%",model$algebraic[[k]]$notes,sep=" "))
    FILETEXT <- paste(FILETEXT,ALGtext,"\n",sep="")
  }
  FILETEXT <- paste(FILETEXT," \n",sep="")

  ###############################################
  # Initial conditions
  ###############################################

  # states
  for (k in seq_along(model$states))
    FILETEXT <- paste(FILETEXT, model$states[[k]]$name, "(0) = ",
                      model$states[[k]]$IC, "\n", sep="")

  # algebraic states
  for (k in seq_along(model$algebraic))
    if (!is.null(model$algebraic[[k]]$name))
      FILETEXT <- paste(FILETEXT, model$algebraic[[k]]$name,
                        "(0) = ", model$algebraic[[k]]$IC, "\n",sep="")

  FILETEXT <- paste(FILETEXT," \n",sep="")

  ###############################################
  # Parameters
  ###############################################

  FILETEXT <- paste(FILETEXT,"********** MODEL PARAMETERS\n\n",sep="")

  for (k in seq_along(model$parameters)) {
    PARtext <- paste(model$parameters[[k]]$name," = ",model$parameters[[k]]$value,sep="")
    type <- model$parameters[[k]]$type
    compartment <- model$parameters[[k]]$compartment
    unittype <- model$parameters[[k]]$unittype
    informationText <- ""

    if (!is.null(type) || !is.null(compartment) || !is.null(unittype)) {
      if (type=="isSpecie" && !is.null(compartment) && (unittype=="amount" || unittype=="concentration"))
        informationText <- paste(" {",type,":",compartment,":",unittype,"}",sep="")
      if (type=="isParameter" && is.null(compartment) && is.null(unittype))
        informationText <- paste(" {",type,"}",sep="")
      if (type=="isCompartment" && is.null(unittype))
        informationText <- paste(" {",type,":",compartment,"}",sep="")
      if (informationText=="") {
        stop(paste("exportTxtAZRmodel: Type information for parameter ",model$parameters[[k]]$name," seems to be wrong."),sep="")
      }
    }

    PARtext <- strtrimM(paste(PARtext,informationText,sep=""))

    if (model$parameters[[k]]$estimate)
      PARtext <- strtrimM(paste(PARtext,"<estimate>",sep=" "))
    if (model$parameters[[k]]$regressor)
      PARtext <- strtrimM(paste(PARtext,"<regressor>",sep=" "))

    if (!is.null(model$parameters[[k]]$notes))
      PARtext <- strtrimM(paste(PARtext,"%",model$parameters[[k]]$notes,sep=" "))
    FILETEXT <- paste(FILETEXT,PARtext,"\n",sep="")
  }
  FILETEXT <- paste(FILETEXT," \n",sep="")

  ###############################################
  # Variables
  ###############################################

  FILETEXT <- paste(FILETEXT,"********** MODEL VARIABLES\n\n",sep="")

  for (k in seq_along(model$variables)) {
    VARtext <- paste(model$variables[[k]]$name," = ",model$variables[[k]]$formula,sep="")
    type <- model$variables[[k]]$type
    compartment <- model$variables[[k]]$compartment
    unittype <- model$variables[[k]]$unittype
    informationText <- ""

    if (!is.null(type) || !is.null(compartment) || !is.null(unittype)) {
      if (type=="isSpecie" && !is.null(compartment) && (unittype=="amount" || unittype=="concentration"))
        informationText <- paste(" {",type,":",compartment,":",unittype,"}",sep="")
      if (type=="isParameter" && is.null(compartment) && is.null(unittype))
        informationText <- paste(" {",type,"}",sep="")
      if (type=="isCompartment" && is.null(unittype))
        informationText <- paste(" {",type,":",compartment,"}",sep="")
      if (informationText=="") {
        stop(paste("exportTxtAZRmodel: Type information for variable ",model$variables[[k]]$name," seems to be wrong."),sep="")
      }
    }

    VARtext <- strtrimM(paste(VARtext,informationText,sep=""))
    if (!is.null(model$variables[[k]]$notes))
      VARtext <- strtrimM(paste(VARtext,"%",model$variables[[k]]$notes,sep=" "))
    FILETEXT <- paste(FILETEXT,VARtext,"\n",sep="")
  }
  FILETEXT <- paste(FILETEXT," \n",sep="")

  ###############################################
  # Reactions
  ###############################################

  FILETEXT <- paste(FILETEXT,"********** MODEL REACTIONS\n\n",sep="")

  for (k in seq_along(model$reactions)) {
    REAtext <- paste(model$reactions[[k]]$name," = ",model$reactions[[k]]$formula,sep="")
    if (model$reactions[[k]]$reversible) REAtext <- paste(REAtext,"{reversible}")
    if (model$reactions[[k]]$fast) REAtext <- paste(REAtext,"{fast}")
    if (!is.null(model$reactions[[k]]$notes))
      REAtext <- strtrimM(paste(REAtext,"%",model$reactions[[k]]$notes,sep=" "))
    FILETEXT <- paste(FILETEXT,REAtext,"\n",sep="")
  }
  FILETEXT <- paste(FILETEXT," \n",sep="")

  ###############################################
  # Functions
  ###############################################

  FILETEXT <- paste(FILETEXT,"********** MODEL FUNCTIONS\n\n",sep="")

  for (k in seq_along(model$functions)) {
    FUNtext <- paste(model$functions[[k]]$name,"(",model$functions[[k]]$arguments,
                     ") = ",model$functions[[k]]$formula,sep="")
    if (!is.null(model$functions[[k]]$notes))
      FUNtext <- strtrimM(paste(FUNtext,"%",model$functions[[k]]$notes,sep=" "))
    FILETEXT <- paste(FILETEXT,FUNtext,"\n",sep="")
  }
  FILETEXT <- paste(FILETEXT," \n",sep="")

  ###############################################
  # Events
  ###############################################

  FILETEXT <- paste(FILETEXT,"********** MODEL EVENTS\n\n",sep="")

  for (k in seq_along(model$events)) {
    EVEtext <- paste(model$events[[k]]$name," = ",model$events[[k]]$trigger,sep="")
    if (getNumberOfEventassignmentsAZRmodel(model,k) > 0) {
      for (k2 in 1:length(model$events[[k]]$assignment))
        EVEtext <- paste(EVEtext,",",model$events[[k]]$assignment[[k2]]$variable,
                         ",",model$events[[k]]$assignment[[k2]]$formula,sep="")
    }
    if (!is.null(model$events[[k]]$notes))
      EVEtext <- strtrimM(paste(EVEtext,"%",model$events[[k]]$notes,sep=" "))
    FILETEXT <- paste(FILETEXT,EVEtext,"\n",sep="")
  }
  FILETEXT <- paste(FILETEXT," ",sep="")

  # Write the file
  readr::write_file(FILETEXT,filename)
}
