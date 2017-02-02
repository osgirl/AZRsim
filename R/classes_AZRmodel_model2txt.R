###############################################################################
###############################################################################
# This file contains the export functions for AZRmodels to TXT  files
# These are not exported and called from AZRexportAZRmodel
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

  if (!is.AZRmodel(model))
    stop("exportTxtAZRmodel: input argument is not an AZRmodel")

  if (is.null(filename))
    filename <- gsub("\\W","",model$name)

  filename <- paste(strrep(filename,".txt",""), ".txt", sep="")

  # Open the file
  fid <- fopen(filename)

  fwrite(fid,paste("********** MODEL NAME\n\n",model$name,"\n",sep=""));

  fwrite(fid,paste("********** MODEL NOTES\n\n",model$notes,"\n",sep=""));

  ###############################################
  # States
  ###############################################

  fwrite(fid,"********** MODEL STATES\n");

  if (getNumberOfStatesAZRmodel(model) > 0) {
    for (k in 1:getNumberOfStatesAZRmodel(model)) {
      ODEtext <- paste("d/dt(",model$states[[k]]$name,") = ",model$states[[k]]$ODE,sep="")
      type <- model$states[[k]]$type;
      compartment <- model$states[[k]]$compartment;
      unittype <- model$states[[k]]$unittype;
      informationText <- ""

      if (!is.null(type) || !is.null(compartment) || !is.null(unittype)) {
        if (type=="isSpecie" && !is.null(compartment) && (unittype=="amount" || unittype=="concentration"))
          informationText <- paste(" {",type,":",compartment,":",unittype,"}",sep="")
        if (type=="isParameter" && is.null(compartment) && is.null(unittype))
          informationText <- paste(" {",type,"}",sep="")
        if (type=="isCompartment" && is.null(unittype))
          informationText <- paste(" {",type,":",compartment,"}",sep="")
        if (informationText=="") {
          fclose(fid)
          stop(paste("exportTxtAZRmodel: Type information for state ",model$states[[k]]$name," seems to be wrong."),sep="")
        }
      }

      if (!is.null(model$states[[k]]$lowConstraint) && !is.null(model$states[[k]]$highConstraint)) {
        constraintsText = paste(" {constraints:[",model$states[[k]]$lowConstraint,",",model$states[[k]]$highConstraint,"]}",sep="")
      } else {
        constraintsText = ""
      }

      ODEtext <- strtrim(paste(ODEtext,informationText,constraintsText,sep=""))
      if (!is.null(model$states[[k]]$notes))
        ODEtext <- strtrim(paste(ODEtext,"%",model$states[[k]]$notes,sep=" "))
      fwrite(fid,ODEtext)
    }
    fwrite(fid," ");
  }

  ###############################################
  # Algebraic Rules
  ###############################################

  if (getNumberOfAlgebraicAZRmodel(model) > 0) {
    for (k in 1:length(model$algebraic)) {
      if (!is.null(model$algebraic[[k]]$name)) {
        ALGtext <- paste("0 = ",model$algebraic[[k]]$formula," : ",model$algebraic[[k]]$name,sep="")
      } else {
        ALGtext <- paste("0 = ",model$algebraic[[k]]$formula,sep="")
      }
      type <- model$algebraic[[k]]$type;
      compartment <- model$algebraic[[k]]$compartment;
      unittype <- model$algebraic[[k]]$unittype;
      informationText <- ""

      if (!is.null(type) || !is.null(compartment) || !is.null(unittype)) {
        if (type=="isSpecie" && !is.null(compartment) && (unittype=="amount" || unittype=="concentration"))
          informationText <- paste(" {",type,":",compartment,":",unittype,"}",sep="")
        if (type=="isParameter" && is.null(compartment) && is.null(unittype))
          informationText <- paste(" {",type,"}",sep="")
        if (type=="isCompartment" && is.null(unittype))
          informationText <- paste(" {",type,":",compartment,"}",sep="")
        if (informationText=="") {
          fclose(fid)
          stop(paste("exportTxtAZRmodel: Type information for algebraic state ",model$algebraic[[k]]$name," seems to be wrong."),sep="")
        }
      }

      ALGtext <- strtrim(paste(ALGtext,informationText,sep=""))
      if (!is.null(model$algebraic[[k]]$notes))
        ALGtext <- strtrim(paste(ALGtext,"%",model$algebraic[[k]]$notes,sep=" "))
      fwrite(fid,ALGtext)
    }
    fwrite(fid," ")
  }

  ###############################################
  # Initial conditions
  ###############################################

  # states
  if (getNumberOfStatesAZRmodel(model) > 0) {
    for (k in 1:length(model$states))
      fwrite(fid,paste(model$states[[k]]$name,"(0) = ",model$states[[k]]$IC,sep=""))
  }

  # algebraic states
  if (getNumberOfAlgebraicAZRmodel(model) > 0) {
    for (k in 1:length(model$algebraic))
      if (!is.null(model$algebraic[[k]]$name))
        fwrite(fid,paste(model$algebraic[[k]]$name,"(0) = ",model$algebraic[[k]]$IC,sep=""))
  }
  fwrite(fid," ")

  ###############################################
  # Parameters
  ###############################################

  fwrite(fid,"********** MODEL PARAMETERS\n");

  if (getNumberOfParametersAZRmodel(model) > 0) {
    for (k in 1:length(model$parameters)) {
      PARtext <- paste(model$parameters[[k]]$name," = ",model$parameters[[k]]$value,sep="")
      type <- model$parameters[[k]]$type;
      compartment <- model$parameters[[k]]$compartment;
      unittype <- model$parameters[[k]]$unittype;
      informationText <- ""

      if (!is.null(type) || !is.null(compartment) || !is.null(unittype)) {
        if (type=="isSpecie" && !is.null(compartment) && (unittype=="amount" || unittype=="concentration"))
          informationText <- paste(" {",type,":",compartment,":",unittype,"}",sep="")
        if (type=="isParameter" && is.null(compartment) && is.null(unittype))
          informationText <- paste(" {",type,"}",sep="")
        if (type=="isCompartment" && is.null(unittype))
          informationText <- paste(" {",type,":",compartment,"}",sep="")
        if (informationText=="") {
          fclose(fid)
          stop(paste("exportTxtAZRmodel: Type information for parameter ",model$parameters[[k]]$name," seems to be wrong."),sep="")
        }
      }

      PARtext <- strtrim(paste(PARtext,informationText,sep=""))

      if (model$parameters[[k]]$estimate)
        PARtext <- strtrim(paste(PARtext,"<estimate>",sep=" "))
      if (model$parameters[[k]]$regressor)
        PARtext <- strtrim(paste(PARtext,"<regressor>",sep=" "))

      if (!is.null(model$parameters[[k]]$notes))
        PARtext <- strtrim(paste(PARtext,"%",model$parameters[[k]]$notes,sep=" "))
      fwrite(fid,PARtext)
    }
  }
  fwrite(fid," ")

  ###############################################
  # Variables
  ###############################################

  fwrite(fid,"********** MODEL VARIABLES\n");

  if (getNumberOfVariablesAZRmodel(model) > 0) {
    for (k in 1:length(model$variables)) {
      VARtext <- paste(model$variables[[k]]$name," = ",model$variables[[k]]$formula,sep="")
      type <- model$variables[[k]]$type;
      compartment <- model$variables[[k]]$compartment;
      unittype <- model$variables[[k]]$unittype;
      informationText <- ""

      if (!is.null(type) || !is.null(compartment) || !is.null(unittype)) {
        if (type=="isSpecie" && !is.null(compartment) && (unittype=="amount" || unittype=="concentration"))
          informationText <- paste(" {",type,":",compartment,":",unittype,"}",sep="")
        if (type=="isParameter" && is.null(compartment) && is.null(unittype))
          informationText <- paste(" {",type,"}",sep="")
        if (type=="isCompartment" && is.null(unittype))
          informationText <- paste(" {",type,":",compartment,"}",sep="")
        if (informationText=="") {
          fclose(fid)
          stop(paste("exportTxtAZRmodel: Type information for variable ",model$variables[[k]]$name," seems to be wrong."),sep="")
        }
      }

      VARtext <- strtrim(paste(VARtext,informationText,sep=""))
      if (!is.null(model$variables[[k]]$notes))
        VARtext <- strtrim(paste(VARtext,"%",model$variables[[k]]$notes,sep=" "))
      fwrite(fid,VARtext)
    }
  }
  fwrite(fid," ")

  ###############################################
  # Reactions
  ###############################################

  fwrite(fid,"********** MODEL REACTIONS\n");

  if (getNumberOfReactionsAZRmodel(model) > 0) {
    for (k in 1:length(model$reactions)) {
      REAtext <- paste(model$reactions[[k]]$name," = ",model$reactions[[k]]$formula,sep="")
      if (model$reactions[[k]]$reversible) REAtext <- paste(REAtext,"{reversible}")
      if (model$reactions[[k]]$fast) REAtext <- paste(REAtext,"{fast}")
      if (!is.null(model$reactions[[k]]$notes))
        REAtext <- strtrim(paste(REAtext,"%",model$reactions[[k]]$notes,sep=" "))
      fwrite(fid,REAtext)
    }
  }
  fwrite(fid," ")

  ###############################################
  # Functions
  ###############################################

  fwrite(fid,"********** MODEL FUNCTIONS\n");

  if (getNumberOfFunctionsAZRmodel(model) > 0) {
    for (k in 1:length(model$functions)) {
      FUNtext <- paste(model$functions[[k]]$name,"(",model$functions[[k]]$arguments,") = ",model$functions[[k]]$formula,sep="")
      if (!is.null(model$functions[[k]]$notes))
        FUNtext <- strtrim(paste(FUNtext,"%",model$functions[[k]]$notes,sep=" "))
      fwrite(fid,FUNtext)
    }
  }
  fwrite(fid," ")

  ###############################################
  # Events
  ###############################################

  fwrite(fid,"********** MODEL EVENTS\n");

  if (getNumberOfEventsAZRmodel(model) > 0) {
    for (k in 1:length(model$events)) {
      EVEtext <- paste(model$events[[k]]$name," = ",model$events[[k]]$trigger,sep="")
      if (getNumberOfEventassignmentsAZRmodel(model,k) > 0) {
        for (k2 in 1:length(model$events[[k]]$assignment))
          EVEtext <- paste(EVEtext,",",model$events[[k]]$assignment[[k2]]$variable,",",model$events[[k]]$assignment[[k2]]$formula,sep="");
      }
      if (!is.null(model$events[[k]]$notes))
        EVEtext <- strtrim(paste(EVEtext,"%",model$events[[k]]$notes,sep=" "))
      fwrite(fid,EVEtext)
    }
  }
  fwrite(fid," ")

  # Close the file
  fclose(fid)
}
