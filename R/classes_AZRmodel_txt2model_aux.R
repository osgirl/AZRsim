###############################################################################
###############################################################################
# This file contains auxiliary functions for the import functions from TXT
# and TXTBC files to AZRmodels
###############################################################################
###############################################################################


###############################################################################
# checkgetNotes - parse out potential comments / notes
###############################################################################
checkgetNotes <- function(textString) {
  # check if comment present
  startNotes <- regexpr("%",textString)
  notesk <- NULL
  if (startNotes[1] != -1) {
    notesk <- strtrim(substr(textString,startNotes[1]+1,nchar(textString)))
    textString <- strtrim(substr(textString,1,startNotes[1]-1))
    if (nchar(notesk)==0) notesk <- NULL
  }
  # create result list
  commentInfo <- list(main=textString,comment=notesk)
  return(commentInfo)
}


###############################################################################
# checkGetFlag - parse out special information such as {reversible}, {fast},
# {estimate}, {regressor}, etc.
###############################################################################
checkGetFlag <- function(textString,flagString) {
  flagPresent <- FALSE
  # check if the "flagString" identifier is present.
  temp <-
  if (!is.null(strlocateall(textString,flagString)$start)) {
    # flagString is present - take it away and
    # set the flag to one, otherwise leave the expression untouched and
    # set it to 0
    textString <- strtrim(strrep(textString,flagString,""))
    flagPresent <- TRUE
  }
  flagInfo <- list(textString=textString,flagPresent=flagPresent)
}


###############################################################################
# checkGetSBMLinfo - parse out potential SBML information
###############################################################################
checkGetSBMLinfo <- function(textString,fctname,elementname) {
  # check if additional information is present ... if yes, cut it out
  infoStart <- strlocateall(textString,"{")
  infoEnd <- strlocateall(textString,"}")
  informationText <- ""
  if ((length(infoStart$start)+length(infoEnd$start))>2) {
    stop(paste(fctname, ": To many curly parentheses in a ",elementname," definition",sep=""))
  }
  if (length(infoStart$start) != length(infoEnd$start)) {
    stop(paste(fctname, ": At least one ",elementname," information not properly defined",sep=""))
  }
  if (!is.null(infoStart$start)) {
    informationText <- strtrim(substr(textString,infoStart$start+1,infoEnd$start-1))
    textString <- strtrim(substr(textString,1,infoStart$start-1))
  }

  type <- NULL
  compartment <- NULL
  unittype <- NULL
  if(nchar(informationText)>0) {
    # explode the information text with ':'
    terms <- strexplode(informationText,':')
    found <- FALSE
    if (tolower(strtrim(terms[1]))=="isparameter") {
      type <- "isParameter"
      found <- TRUE
    }
    if (tolower(strtrim(terms[1]))=="iscompartment") {
      type <- "isCompartment"
      if (length(terms)==1) terms[2] = ""
      compartment <- strtrim(terms[2])
      found <- TRUE
    }
    if (tolower(strtrim(terms[1]))=="isspecie") {
      type <- "isSpecie"
      if (length(terms)!=3)
        stop(paste(fctname, ": Error in ",elementname," isSpecie SBML information",sep=""))
      compartment <- strtrim(terms[2])
      unittype <- strtrim(terms[3])
      found <- TRUE
    }
    if (!found)
      stop(paste(fctname, ": Error in ",elementname," SBML information",sep=""))
  }
  SBMLinfo <- list(textString=textString,type=type,compartment=compartment,unittype=unittype)
  return(SBMLinfo)
}


###############################################################################
# getParameters
###############################################################################
getParameters <- function(model,model_parameters) {

  # run through the variables and process them
  if (!is.null(model_parameters)) {
    for (k in 1:length(model_parameters)) {
      parameterString <- strtrim(model_parameters[k])

      # Parse comments / notes
      commentInfo     <- checkgetNotes(parameterString)
      parameterString <- commentInfo$main
      notesk          <- commentInfo$comment

      # check if the "<estimate>" identifier is present.
      flagInfo        <- checkGetFlag(parameterString,"<estimate>")
      parameterString <- flagInfo$textString
      estimateFlag    <- flagInfo$flagPresent

      # check if the "<regressor>" identifier is present.
      flagInfo        <- checkGetFlag(parameterString,"<regressor>")
      parameterString <- flagInfo$textString
      regressorFlag   <- flagInfo$flagPresent

      # Check if both estimate and regressor flag given
      if (estimateFlag && regressorFlag)
        stop("getParameters: parameter with both <estimate> and <regressor> flag")

      # Parse SBML related information
      SBMLinfo        <- checkGetSBMLinfo(parameterString,"getParameters","parameter")
      typek           <- SBMLinfo$type
      compartmentk    <- SBMLinfo$compartment
      unittypek       <- SBMLinfo$unittype
      parameterString <- SBMLinfo$textString

      # extract the parameter name
      temp <- regexpr("=", parameterString)
      test <- strtrim(substr(parameterString,1,(temp[1]-1)))
      # check if parameter name given
      if (nchar(test) == 0) {
        stop("getParameters: At least one parameter name not given.")
      }
      namek <- strremWhite(test)

      # extract the parameter value
      valuek = strtrim(substr(parameterString,(temp+1),nchar(parameterString)))

      # check if parameter value given
      if (nchar(valuek) == 0) {
        stop("getParameters: At least one parameter definition not given.")
      }

      # add parameter to model
      model <- addParameterAZRmodel(model,name=namek,value=as.numeric(valuek),
                                    notes=notesk,type=typek,compartment=compartmentk,
                                    unittype=unittypek,estimate=estimateFlag,regressor=regressorFlag)
    }
  }
  return(model)
}


###############################################################################
# getVariables
###############################################################################
getVariables <- function(model,model_variables) {

  # run through the variables and process them
  if (!is.null(model_variables)) {
    for (k in 1:length(model_variables)) {
      variableString <- strtrim(model_variables[k])

      # Parse comments / notes
      commentInfo    <- checkgetNotes(variableString)
      variableString <- commentInfo$main
      notesk         <- commentInfo$comment

      # Parse SBML related information
      SBMLinfo        <- checkGetSBMLinfo(variableString,"getVariables","variable")
      typek           <- SBMLinfo$type
      compartmentk    <- SBMLinfo$compartment
      unittypek       <- SBMLinfo$unittype
      variableString  <- SBMLinfo$textString

      # extract the variable name
      temp <- regexpr("=", variableString)
      test <- strtrim(substr(variableString,1,(temp[1]-1)))
      # check if variable name given
      if (nchar(test) == 0) {
        stop("getVariables: At least one variable name not given.")
      }
      namek <- strremWhite(test)

      # extract the variable expression
      formulak = strtrim(substr(variableString,(temp+1),nchar(variableString)))

      # check if variable expression given
      if (nchar(formulak) == 0) {
        stop("getVariables: At least one variable definition not given.")
      }

      # add variable to model
      model <- addVariableAZRmodel(model,name=namek,formula=formulak,notes=notesk,type=typek,compartment=compartmentk,unittype=unittypek)
    }
  }
  return(model)
}

###############################################################################
# getFunctions
###############################################################################
getFunctions <- function(model,model_functions) {

  if (!is.null(model_functions)) {
    # run through the functions and process them
    for (k in 1:length(model_functions)) {

      functionString = strtrim(model_functions[k])

      # Parse comments / notes
      commentInfo    <- checkgetNotes(functionString)
      functionString <- commentInfo$main
      notesk         <- commentInfo$comment

      # function name
      temp <- regexpr("\\(", functionString)
      namek <- strtrim(substr(functionString,1,(temp[1]-1)))

      # check if function name given
      if (nchar(namek) == 0) {
        stop("getFunctions: At least one function name not given.")
      }
      namek <- strremWhite(namek)

      # function arguments
      temp2 <- regexpr("\\)", functionString)
      test <- strtrim(substr(functionString,(temp[1]+1),(temp2[1]-1)))
      # check if function arguments are given
      if (nchar(test) == 0) {
        stop("getFunctions: At least for one function no arguments given.")
      }
      argumentsk <- strremWhite(test)

      # extract the formula
      temp3 <- regexpr("=", functionString)
      formulak <- strtrim(substr(functionString,(temp3[1]+1),nchar(functionString)))

      # check if function formula given
      if(nchar(formulak) == 0) {
        stop("getFunctions: At least for one function no formula given.")
      }

      # add info about the function parts for kth function into the function lists
      model <- addFunctionAZRmodel(model,name=namek,arguments=argumentsk,formula=formulak,notes=notesk)
    }
  }
  return(model)
}

###############################################################################
# getEvents
###############################################################################
getEvents <- function(model,model_events) {

  # run through the events and process them
  if (!is.null(model_events)) {
    for (k in 1:length(model_events)) {
      eventString <- strtrim(model_events[k])

      # Parse comments / notes
      commentInfo <- checkgetNotes(eventString)
      eventString <- commentInfo$main
      notesk      <- commentInfo$comment

      # extract the event name
      temp <- regexpr("=", eventString)
      namek <- strtrim(substr(eventString,1,(temp[1]-1)))
      # check if event name given
      if (nchar(namek) == 0) {
        stop("getEvents: At least one event name not given.")
      }
      # get the right hand side
      eventRHS <- strtrim(substr(eventString,temp[1]+1,nchar(eventString)))
      # decompose the eventRHS into its comma separated elements
      # taking into account parentheses
      elementsRHS <- strexplodePC(eventRHS)
      # check number of elements
      if ((length(elementsRHS) < 3) | ((length(elementsRHS) %% 2) == 0)) {
        stop("getEvents: At least one event has no full information given.")
      }
      # first element is assumed to be the trigger function
      triggerk <- strremWhite(elementsRHS[1])
      # Add event to the model
      model <- addEventAZRmodel(model,name=namek,trigger=strremWhite(triggerk),notes=notesk)
      # Add event assignments
      for (k2 in seq(2,length(elementsRHS),2)) {
        model <- addEventAssignmentAZRmodel(model,eventindex=k,variable=strremWhite(elementsRHS[k2]),formula=strremWhite(elementsRHS[k2+1]))
      }
    }
  }
  return(model)
}
