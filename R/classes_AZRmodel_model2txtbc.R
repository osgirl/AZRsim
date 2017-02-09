###############################################################################
###############################################################################
# This file contains the export functions for AZRmodels to TXTBC files
# These are not exported and called from AZRexportAZRmodel
###############################################################################
###############################################################################

###############################################################################
# exportTxtBcAZRmodel: exports AZRmodel as .txtbc file
###############################################################################
# Export of AZRmodel to .txtbc file
#
# Export of an AZRmodel to a .txtbc file. Biochemical reaction notation is used
# where possible, otherwise ODEs.
#
# @param model An AZRmodel to be exported.
# @param filename Full path with filename to export the model to (.txtbc will be added)
# @return None

exportTxtBcAZRmodel <- function (model, filename=NULL) {

  if (!is.AZRmodel(model))
    stop("exportTxtBcAZRmodel: input argument is not an AZRmodel")

  if (is.null(filename))
    filename <- gsub("\\W","",model$name)

  filename <- paste(strrepM(filename,".txtbc",""), ".txtbc", sep="")

  # Initialize the FILETEXT
  FILETEXT <- ""

  FILETEXT <- paste(FILETEXT,"********** MODEL NAME\n\n",model$name,"\n\n",sep="")

  FILETEXT <- paste(FILETEXT,"********** MODEL NOTES\n\n",model$notes,"\n\n",sep="")

  ###############################################
  # Get information about the stoichiometry
  ###############################################
  # first check if the complete stoichiometric matrix can be determined.
  # otherwise a biochemical representation is not fully possible and some
  # states still need to be defined by differential equations.
  # determine the names of that states that need to be described by ODEs
  stoichInfo <- stoichiometryAZRmodel(model,raw=FALSE)
  N <- stoichInfo$N
  stateNamesAll <- getAllStatesAZRmodel(model)$statenames
  stateNamesBC <- stoichInfo$statenames
  stateNamesODE <- setdiff(stateNamesAll,stateNamesBC)

  ###############################################
  # PROCESS THE states
  ###############################################

  FILETEXT <- paste(FILETEXT,"********** MODEL STATE INFORMATION\n\n",sep="")

  # first the definition of the states+ODEs that need to be described by ODEs
  if (length(stateNamesODE) > 0) {
    for (k in 1:length(stateNamesODE)) {
      index <- strmatch(stateNamesAll,stateNamesODE[k])
      ODEtext <- paste("d/dt(",model$states[[index]]$name,") = ",model$states[[index]]$ODE,sep="")
      FILETEXT <- paste(FILETEXT,ODEtext,"\n",sep="")
    }
    FILETEXT <- paste(FILETEXT," \n",sep="")
  }

  # WRITE OUT ALGEBRAIC RULES
  if (getNumberOfAlgebraicAZRmodel(model) > 0) {
    for (k in 1:length(model$algebraic)) {
      if (!is.null(model$algebraic[[k]]$name)) {
        ALGtext <- paste("0 = ",model$algebraic[[k]]$formula," : ",model$algebraic[[k]]$name,sep="")
      } else {
        ALGtext <- paste("0 = ",model$algebraic[[k]]$formula,sep="")
      }
      FILETEXT <- paste(FILETEXT,ALGtext,"\n",sep="")
    }
  }
  FILETEXT <- paste(FILETEXT," \n",sep="")

  # WRITE OUT INITIAL CONDITIONS
  # definition of the initial conditions of the states.
  # additionally for each component optional additional information
  # is processed and written behind the initial conditions.
  # states.type
  # states.compartment
  # states.unittype

  # now construct the states IC text
  # states
  if (getNumberOfStatesAZRmodel(model) > 0) {
    for (k in 1:getNumberOfStatesAZRmodel(model)) {
      ICtext <- paste(model$states[[k]]$name,"(0) = ",model$states[[k]]$IC,sep="")

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
          fclose(fid)
          stop(paste("exportTxtBcAZRmodel: Type information for state ",model$states[[k]]$name," seems to be wrong."),sep="")
        }
      }

      if (!is.null(model$states[[k]]$lowConstraint) && !is.null(model$states[[k]]$highConstraint)) {
        constraintsText = paste(" {constraints:[",model$states[[k]]$lowConstraint,",",model$states[[k]]$highConstraint,"]}",sep="")
      } else {
        constraintsText = ""
      }

      ICtext <- strtrimM(paste(ICtext,informationText,constraintsText,sep=""))

      if (!is.null(model$states[[k]]$notes))
        ICtext <- strtrimM(paste(ICtext,"%",model$states[[k]]$notes,sep=" "))

      FILETEXT <- paste(FILETEXT,ICtext,"\n",sep="")

    }
    FILETEXT <- paste(FILETEXT," \n",sep="")
  }

  # do the same for algebraic variable initial conditions
  if (getNumberOfAlgebraicAZRmodel(model) > 0) {
    for (k in 1:length(model$algebraic)) {
      if (!is.null(model$algebraic[[k]]$name)) {
        ICtext <- paste(model$algebraic[[k]]$name,"(0) = ",model$algebraic[[k]]$IC,sep="")
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
            fclose(fid)
            stop(paste("exportTxtBcAZRmodel: Type information for algebraic state ",model$algebraic[[k]]$name," seems to be wrong."),sep="")
          }
        }

        ICtext <- strtrimM(paste(ICtext,informationText,sep=""))
        if (!is.null(model$algebraic[[k]]$notes))
          ICtext <- strtrimM(paste(ICtext,"%",model$algebraic[[k]]$notes,sep=" "))
        FILETEXT <- paste(FILETEXT,ICtext,"\n",sep="")
      }
    }
    FILETEXT <- paste(FILETEXT," \n",sep="")
  }

  ###############################################
  # Parameters
  ###############################################

  FILETEXT <- paste(FILETEXT,"********** MODEL PARAMETERS\n\n",sep="")

  if (getNumberOfParametersAZRmodel(model) > 0) {
    for (k in 1:length(model$parameters)) {
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
          fclose(fid)
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
  }
  FILETEXT <- paste(FILETEXT," \n",sep="")

  ###############################################
  # Variables
  ###############################################

  FILETEXT <- paste(FILETEXT,"********** MODEL VARIABLES\n\n",sep="")

  if (getNumberOfVariablesAZRmodel(model) > 0) {
    for (k in 1:length(model$variables)) {
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
          fclose(fid)
          stop(paste("exportTxtAZRmodel: Type information for variable ",model$variables[[k]]$name," seems to be wrong."),sep="")
        }
      }

      VARtext <- strtrimM(paste(VARtext,informationText,sep=""))
      if (!is.null(model$variables[[k]]$notes))
        VARtext <- strtrimM(paste(VARtext,"%",model$variables[[k]]$notes,sep=" "))
      FILETEXT <- paste(FILETEXT,VARtext,"\n",sep="")

    }
  }
  FILETEXT <- paste(FILETEXT," \n",sep="")

  ###############################################
  # PROCESS THE REACTIONS
  ###############################################

  FILETEXT <- paste(FILETEXT,"********** MODEL REACTIONS\n\n",sep="")

  # Use the stoichiometric matrix determined above (by eventually adding
  # species that are defined by variables)
  # get reaction names, kinetics, and reversibility flag
  reacInfo <- getAllReactionsAZRmodel(model)

  # cycle through the columns of the stoichiometric matrix N to build the
  # reaction expressions
  if (length(N)!=0) {
    for (k1 in 1:dim(N)[2]) {
      Ncol <- N[,k1]
      # first get the substrates by finding the negative elements in Ncol
      substrateIndices <- unname(which(Ncol < 0))
      # then get the products by finding the positive elements in Ncol
      productIndices <- unname(which(Ncol > 0))
      # determine the needed information
      reactionName <- reacInfo$reacnames[k1]
      reactionRevFlag <- reacInfo$reacreversible[k1]
      reactionFastFlag <- reacInfo$reacfast[k1]
      reactionFormula <- reacInfo$reacformulas[k1]
      substrateNames <- stateNamesBC[substrateIndices]
      productNames <- stateNamesBC[productIndices]
      substrateStoichiometries <- abs(Ncol[substrateIndices])
      productStoichiometries <- abs(Ncol[productIndices])
      # if reversible split up the reaction rate in two parts. if this is not
      # possible, issue a warning and set the reaction as irreversible
      if (reactionRevFlag) {
        irreversibleRates <- strexplodePC(reactionFormula,'-')
        if (length(irreversibleRates) != 2) {
          # Need to have two parts that are separated by a '-' sign. (Forward
          # first, then reverse reaction kinetics).
          reactionRevFlag <- FALSE
        } else {
          reactionForward <- irreversibleRates[1]
          reactionReverse <- irreversibleRates[2]
        }
      }
      # format the output of the reaction text
      reacText <- ""
      # first the reaction expression, e.g. 2*A + 4*C => 3*B
      # the substrates
      if (length(substrateNames) > 0) {
        if (substrateStoichiometries[1] != 1) {
          reacText <- paste(reacText,substrateStoichiometries[1],"*",substrateNames[1],sep="")
        } else {
          reacText <- paste(reacText,substrateNames[1],sep="")
        }
      }
      if (length(substrateNames) > 1) {
        for (k2 in 2:length(substrateNames)) {
          if (substrateStoichiometries[k2] != 1) {
            reacText <- paste(reacText,"+",substrateStoichiometries[k2],"*",substrateNames[k2],sep="")
          } else {
            reacText <- paste(reacText,"+",substrateNames[k2],sep="")
          }
        }
      }
      # the reaction equation sign
      if (reactionRevFlag) {
        reacText <- paste(reacText," <=> ",sep="")
      } else {
        reacText <- paste(reacText," => ",sep="")
      }
      # the products
      if (length(productNames) > 0) {
        if (productStoichiometries[1] != 1) {
          reacText <- paste(reacText,productStoichiometries[1],"*",productNames[1],sep="")
        } else {
          reacText <- paste(reacText,productNames[1],sep="")
        }
      }
      if (length(productNames) > 1) {
        for (k2 in 2:length(productNames)) {
          if (productStoichiometries[k2] != 1) {
            reacText <- paste(reacText,"+",productStoichiometries[k2],"*",productNames[k2],sep="")
          } else {
            reacText <- paste(reacText,"+",productNames[k2],sep="")
          }
        }
      }
      # separator and reaction name
      reacText <- paste(reacText," : ",reactionName,sep="")
      # fast flag
      if (reactionFastFlag) reacText <- paste(reacText," {fast}",sep="")
      # notes
      if (!is.null(model$reactions[[k1]]$notes)) reacText <- paste(reacText," % ",model$reactions[[k1]]$notes,sep="")
      # new line
      reacText <- paste(reacText,"\n",sep="")
      # now the reaction rate expression(s)
      if (!reactionRevFlag) {
        reacText <- paste(reacText,"\tvf = ",reactionFormula,"\n",sep="")
      } else {
        reacText <- paste(reacText,"\tvf = ",reactionForward,"\n\tvr = ",reactionReverse,"\n",sep="")
      }
      # write out reaction text
      FILETEXT <- paste(FILETEXT,reacText,"\n",sep="")

    }
  } else {
    FILETEXT <- paste(FILETEXT," \n",sep="")
  }

  ###############################################
  # Functions
  ###############################################

  FILETEXT <- paste(FILETEXT,"********** MODEL FUNCTIONS\n\n",sep="")

  if (getNumberOfFunctionsAZRmodel(model) > 0) {
    for (k in 1:length(model$functions)) {
      FUNtext <- paste(model$functions[[k]]$name,"(",model$functions[[k]]$arguments,") = ",model$functions[[k]]$formula,sep="")
      if (!is.null(model$functions[[k]]$notes))
        FUNtext <- strtrimM(paste(FUNtext,"%",model$functions[[k]]$notes,sep=" "))
      FILETEXT <- paste(FILETEXT,FUNtext,"\n",sep="")

    }
  }
  FILETEXT <- paste(FILETEXT," \n",sep="")

  ###############################################
  # Events
  ###############################################

  FILETEXT <- paste(FILETEXT,"********** MODEL EVENTS\n\n",sep="")

  if (getNumberOfEventsAZRmodel(model) > 0) {
    for (k in 1:length(model$events)) {
      EVEtext <- paste(model$events[[k]]$name," = ",model$events[[k]]$trigger,sep="")
      if (getNumberOfEventassignmentsAZRmodel(model,k) > 0) {
        for (k2 in 1:length(model$events[[k]]$assignment))
          EVEtext <- paste(EVEtext,",",model$events[[k]]$assignment[[k2]]$variable,",",model$events[[k]]$assignment[[k2]]$formula,sep="")
      }
      if (!is.null(model$events[[k]]$notes))
        EVEtext <- strtrimM(paste(EVEtext,"%",model$events[[k]]$notes,sep=" "))
      FILETEXT <- paste(FILETEXT,EVEtext,"\n",sep="")

    }
  }
  FILETEXT <- paste(FILETEXT," ",sep="")

  # Write the file
  filewrite(FILETEXT,filename)
}
