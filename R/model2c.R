###############################################################################
###############################################################################
# This file contains the export functions for AZRmodels to C files
# to be used with the CVODES integrator
###############################################################################
###############################################################################

###############################################################################
# exportCcodeAZRmodel: exports AZRmodel as .c file
###############################################################################
# Export of AZRmodel as .c file
#
# Export of an AZRmodel to C code, allowing to simulate the model with the
# Sundials CVODES integrator.
#
# @param model An AZRmodel to be exported to C code
# @param filename Full path with filename to export the model to (.c will be added)
# @return None

exportCcodeAZRmodel <- function (model, filename=NULL) {

  ###############################################
  # Check input arguments
  ###############################################
  if (!AZRsim::is_azrmod(model))
    stop("exportCcodeAZRmodel: input argument is not an AZRmodel")

  if (is.null(filename))
    filename <- tempfile()

  #!!!!!!!!!!!!!!!! NEED TO FIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!! NEED TO FIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # ###############################################
  # # ADD PIECEWISE TRIGGERS AS EVENTS
  # ###############################################
  # model <- addpiecewiseeventsIQM(model)
  #!!!!!!!!!!!!!!!! NEED TO FIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!! NEED TO FIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ###############################################
  # Generate .h and .c filenames
  ###############################################
  fileInfo  <- fileparts(filename)
  filenameC <- paste(paste(fileInfo$pathname,fileInfo$filename,sep="/"),"c",sep=".")

  ###############################################
  # GET THE MODELS ELEMENTS
  ###############################################
  stateInfo <- get_all_states(model)
  paramInfo <- get_all_parameters(model)
  varInfo   <- get_all_variables(model)
  reacInfo  <- get_all_reactions(model)
  funcInfo  <- get_all_functions(model)
  eveInfo   <- get_all_events(model)

  ###############################################
  # GET NUMBER OF THE MODELS ELEMENTS
  ###############################################
  NRSTATES      <- len_states(model)
  NRPARAMETERS  <- len_parameters(model)
  NRVARIABLES   <- len_variables(model)
  NRREACTIONS   <- len_reactions(model)
  NRFUNCTIONS   <- len_functions(model)
  NREVENTS      <- len_events(model)

  ###############################################
  # DEAL WITH THE FORMULAS
  ###############################################
  # For C the double representation needs to be made (1->1.0) ...
  # Furthermore, several other things need to be fixed.
  stateInfo$stateODEs   <- dealFormulas(stateInfo$stateODEs)
  varInfo$varformulas   <- dealFormulas(varInfo$varformulas)
  reacInfo$reacformulas <- dealFormulas(reacInfo$reacformulas)
  funcInfo$funcformulas <- dealFormulas(funcInfo$funcformulas)
  eveInfo$evetriggers   <- dealFormulas(eveInfo$evetriggers)
  if (length(eveInfo$eveformulas) > 0) {
    for (k in 1:length(eveInfo$eveformulas)) {
      eveInfo$eveformulas[[k]] <- dealFormulas(eveInfo$eveformulas[[k]])
    }
  }

  ###############################################
  # WRITE THE MODEL C FILE
  ###############################################
  # Initialize the FILETEXT
  FILETEXT <- ""

  FILETEXT <- paste(FILETEXT,'#include <stddef.h>\n',sep="")
  FILETEXT <- paste(FILETEXT,'#include <stdarg.h>\n',sep="")
  FILETEXT <- paste(FILETEXT,'#include <math.h>\n\n',sep="")

  FILETEXT <- paste(FILETEXT,'#include "cvodes_i.h"\n',sep="")
  FILETEXT <- paste(FILETEXT,'#include "model/splineaddon.h"\n',sep="")
  FILETEXT <- paste(FILETEXT,'#include "model/mathaddon.h"\n',sep="")
  FILETEXT <- paste(FILETEXT,'#include "model/kineticformulas.h"\n',sep="")
  FILETEXT <- paste(FILETEXT,' \n',sep="")

  # First define the functions
  if (NRFUNCTIONS > 0) {
    for (k in 1:NRFUNCTIONS) {
      # Start declaration
      textPiece <- paste("static double ",funcInfo$funcnames[k],"(",sep="")
      # Write arguments
      arguments <- strexplode(funcInfo$funcarguments[k])
      for (k2 in 1:length(arguments)) {
        if (k2 < length(arguments)) {
          textPiece <- paste(textPiece,"double ",arguments[k2],",",sep="")
        } else {
          textPiece <- paste(textPiece,"double ",arguments[k2],")",sep="")
        }
      }
      # Write out
      FILETEXT <- paste(FILETEXT,textPiece,'\n',sep="")
      # write forumla and return
      FILETEXT <- paste(FILETEXT,'{\n',sep="")
      FILETEXT <- paste(FILETEXT,paste("    return ",funcInfo$funcformulas[k],";",sep=""),'\n',sep="")
      FILETEXT <- paste(FILETEXT,'}\n\n',sep="")
    }
  }

  # Define the model function
  FILETEXT <- paste(FILETEXT,'void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)\n',sep="")
  FILETEXT <- paste(FILETEXT,'{\n',sep="")

  FILETEXT <- outputDeclarationData(FILETEXT,stateInfo$statenames)
  FILETEXT <- outputDeclarationData(FILETEXT,paramInfo$paramnames)
  FILETEXT <- outputDeclarationData(FILETEXT,varInfo$varnames)
  FILETEXT <- outputDeclarationData(FILETEXT,reacInfo$reacnames)

  if (NREVENTS > 0) {
    eventassignn <- c()
    for (k in 1:NREVENTS) {
      for (k2 in 1:length(eveInfo$eveformulas[[k]])) {
        eventassignn <- c(eventassignn,paste('eventassign_',k,'_',k2,sep=""))
      }
    }
    FILETEXT <- outputDeclarationData(FILETEXT,eventassignn)
  }

  FILETEXT <- paste(FILETEXT,' \n',sep="")

  if (NRSTATES>0) {
    for (k in 1:NRSTATES) {
      FILETEXT <- paste(FILETEXT,paste('    ',stateInfo$statenames[k],' = stateVector[',k-1,'];',sep=""),'\n',sep="")
    }
    FILETEXT <- paste(FILETEXT,' \n',sep="")
  }

  if (NRPARAMETERS>0) {
    for (k in 1:NRPARAMETERS) {
      FILETEXT <- paste(FILETEXT,paste('    ',paramInfo$paramnames[k],' = paramdataPtr->parametervector[',k-1,']; /* ',paramInfo$paramvalues[k],' */',sep=""),'\n',sep="")
    }
    FILETEXT <- paste(FILETEXT,' \n',sep="")
  }

  if (NRVARIABLES>0) {
    for (k in 1:NRVARIABLES) {
      FILETEXT <- paste(FILETEXT,paste('    ',varInfo$varnames[k],' = ',varInfo$varformulas[k],';',sep=""),'\n',sep="")
    }
    FILETEXT <- paste(FILETEXT,' \n',sep="")
  }

  if (NRREACTIONS>0) {
    for (k in 1:NRREACTIONS) {
      FILETEXT <- paste(FILETEXT,paste('    ',reacInfo$reacnames[k],' = ',reacInfo$reacformulas[k],';',sep=""),'\n',sep="")
    }
    FILETEXT <- paste(FILETEXT,' \n',sep="")
  }

  if (NREVENTS>0) {
    for (k in 1:NREVENTS) {
      for (k2 in 1:length(eveInfo$evevariables[[k]])) {
        namevar <- paste('eventassign_',k,'_',k2,sep="")
        FILETEXT <- paste(FILETEXT,paste('    ',namevar,' = ',eveInfo$eveformulas[[k]][k2],';',sep=""),'\n',sep="")
      }
    }
    FILETEXT <- paste(FILETEXT,' \n',sep="")
  }

  FILETEXT <- paste(FILETEXT,'    if (DOflag == DOFLAG_DDT) {\n',sep="")
  if (NRSTATES>0) {
    for (k in 1:NRSTATES) {
      FILETEXT <- paste(FILETEXT,paste('        DDTvector[',k-1,'] = ',stateInfo$stateODEs[k],';',sep=""),'\n',sep="")
    }
  }
  FILETEXT <- paste(FILETEXT,'    } else if (DOflag == DOFLAG_VARREAC) {\n',sep="")
  if (NRVARIABLES>0) {
    for (k in 1:NRVARIABLES) {
      FILETEXT <- paste(FILETEXT,paste('        variableVector[',k-1,'] = ',varInfo$varnames[k],';',sep=""),'\n',sep="")
    }
  }
  if (NRREACTIONS>0) {
    for (k in 1:NRREACTIONS) {
      FILETEXT <- paste(FILETEXT,paste('        reactionVector[',k-1,'] = ',reacInfo$reacnames[k],';',sep=""),'\n',sep="")
    }
  }
  FILETEXT <- paste(FILETEXT,'    } else if (DOflag == DOFLAG_EVENTS) {\n',sep="")
  if (NREVENTS>0) {
    for (k in 1:NREVENTS) {
      tExpr <- getTriggerExpression(eveInfo$evetriggers[k])
      FILETEXT <- paste(FILETEXT,paste('        gout[',k-1,'] = ',tExpr,';',sep=""),'\n',sep="")
    }
  }
  FILETEXT <- paste(FILETEXT,'    } else if (DOflag == DOFLAG_EVENTASSIGN) {\n',sep="")

  if (NREVENTS>0) {
    for (k in 1:NREVENTS) {
      FILETEXT <- paste(FILETEXT,paste('        if (eventVector[',k-1,'] == 1 && gout[',k-1,'] < 0) {',sep=""),'\n',sep="")
      FILETEXT <- paste(FILETEXT,'            DDTvector[0] = 1;\n',sep="")
      vars <- eveInfo$evevariables[[k]]
      for (k2 in 1:length(vars)) {
        index <- strmatch(stateInfo$statenames,vars[k2])-1
        if (length(index) > 0) {
          FILETEXT <- paste(FILETEXT,paste('            stateVector[',index,'] = eventassign_',k,'_',k2,';',sep=""),'\n',sep="")
        }
        index <- strmatch(paramInfo$paramnames,vars[k2])-1
        if (length(index) > 0) {
          FILETEXT <- paste(FILETEXT,paste('            paramdataPtr->parametervector[',index,'] = eventassign_',k,'_',k2,';',sep=""),'\n',sep="")
        }
      }
      FILETEXT <- paste(FILETEXT,'        }\n',sep="")
    }
  }
  FILETEXT <- paste(FILETEXT,'    }\n',sep="")
  FILETEXT <- paste(FILETEXT,'}\n',sep="")
  FILETEXT <- paste(FILETEXT,' ');

  # Write the file
  readr::write_file(FILETEXT,filenameC)
}

##########################################################################
# DEAL WITH THE FORMULAS
# For C the double representation needs to be made (1->1.0) ...
# Furthermore, several other things need to be fixed.
##########################################################################
dealFormulas <- function (formulaArray) {

  if (length(formulaArray)==0) {
    return(formulaArray)
  }

  # Replace names by adding "AZR" at the end
  oldElements = c('\\bnthroot\\b','\\band\\b','\\bor\\b','\\babs\\b','\\bindexmax\\b','\\bmin\\b','\\bmax\\b','\\bpiecewise\\b','\\binterpcs\\b')
  newElements = c('nthrootAZR','andAZR','orAZR','absAZR','indexmaxAZR','minAZR','maxAZR','piecewiseAZR','interpcsAZR')
  for (k in 1:length(oldElements)) formulaArray <- gsub(oldElements[k],newElements[k],formulaArray)

  # Replace name of power operator (power(->pow()
  formulaArray <- gsub('\\bpower\\(','pow(',formulaArray)

  # Handle additional things for each formula in array
  for (k in 1:length(formulaArray)) {

    formula <- formulaArray[k]

    # handle interp0 -> piecewise
    formula <- exchangeInterp0(formula)
    # handle interp1 -> piecewise
    formula <- exchangeInterp1(formula)

    # handle interpcs: changing the syntax
    formula <- exchangeInterpcs(formula)

    # handle the power operator
    formula <- convertPowerOperator(formula)

    # fix the c notation of doubles
    formula <- gsub("(\\b[0-9.]+)","\\1.0",formula)
    formula <- gsub("\\.0\\.",".",formula)
    formula <- gsub("(\\.[0-9]+)\\.0","\\1",formula)
    formula <- gsub("(E\\.0)","E",formula)
    formula <- gsub("([0-9]E-\\d*)\\.0","\\1",formula)
    formula <- gsub("([0-9]E\\+\\d*)\\.0","\\1",formula)
    formula <- gsub("([0-9]E\\d*)\\.0","\\1",formula)

    # Add number of variable input arguments to function calls as first
    # input argument (for the functions, defined above)
    checkElements = c('\\<indexmaxAZR\\>','\\<minAZR\\>','\\<maxAZR\\>','\\<andAZR\\>','\\<orAZR\\>','\\<piecewiseAZR\\>','\\<interpcsAZR\\>')
    for (k1 in 1:length(checkElements)) {
      index <- unlist(gregexpr(pattern=checkElements[k1],formula))
      if (length(index) != 0) {
        for (k2 in 1:length(index)) {
          if (index[k2] != -1) {
            indexStart <- index[k2]+nchar(checkElements[k1])-4
            indexEnd <- indexStart
            parOpen <- 1
            while (parOpen != 0) {
              indexEnd <- indexEnd + 1
              if (substr(formula,indexEnd,indexEnd) == '(') {
                parOpen <- parOpen + 1
              } else {
                if (substr(formula,indexEnd,indexEnd) == ')') {
                  parOpen <- parOpen - 1
                }
              }
            }
            command <- checkElements[k1]
            oldarguments <- substr(formula,indexStart+1,indexEnd-1)
            oldargumentsReplace <- oldarguments
            # For all C interpolation functions the first input argument needs to be
            # reverted to integer.
            if (command %in% c("\\<interpcsAZR\\>")) {
              oldarguments <- sub("([0-9]+).0","\\1",oldarguments)
            }
            newargstring <- paste(length(strexplodePC(oldarguments)),",",oldarguments,sep="")
            oldrep  <- paste(substr(command,3,nchar(command)-2), '(', oldargumentsReplace, ')', sep="")
            newrep <- paste(substr(command,3,nchar(command)-2), '(', newargstring, ')', sep="")
            formula <- strrepM(formula,oldrep,newrep)
            index <- index + nchar(newargstring)-nchar(oldarguments)
          }
        }
      }
    }

    formulaArray[k] <- formula
  }

  return(formulaArray)
}

###############################################
# OUTPUT DECLARATION DATA
###############################################
outputDeclarationData <- function(FILETEXT,data) {
  NR <- length(data)
  if (NR==0) return(FILETEXT)

  # Construct
  textPiece <- ""
  nrperrow <- 0
  for (k in 1:NR) {
    if (nrperrow==0) {
      textPiece <- paste(textPiece,'    double ',sep="")
    }
    if (k<NR && nrperrow<20-1) {
      textPiece <- paste(textPiece,data[k],",",sep="")
    } else {
      textPiece <- paste(textPiece,data[k],";",sep="")
    }
    nrperrow <- nrperrow + 1
    if (nrperrow == 20) {
      textPiece <- paste(textPiece,"\n",sep="")
      nrperrow <- 0
    }
  }
  # Return
  FILETEXT <- paste(FILETEXT,textPiece,"\n",sep="")
  return(FILETEXT)
}

###############################################
# DETERMINE THE TRIGGER EXPRESSION
###############################################
getTriggerExpression <- function (data) {
  tExpr <- paste(data,"-0.5",sep="")
  return(tExpr)
}

##########################################################################
# HANDLE INTERP0 (LOOKUP TABLE W/ ZERO-ORDER INTERPOLATION)
##########################################################################
exchangeInterp0 <- function (text) {
  # Check if interp0 present
  if (is.null(strlocateall(text,"interp0(")$start)) {
    return(text)
  }
  # Check if present multiple times
  if (length(strlocateall(text,"interp0(")$start) > 1) {
    stop('exchangeInterp0: the "interp0" function is only allowed to be present once in each formula.')
  }

  # Present once => handle with piecewise
  textnew <- text

  # Get starting index of things inside
  indexstart <- strlocateall(text,"interp0(")$end+1

  # Cut out the content using the parentheses
  # Get the end index of the k-th statement (closing parentheses
  # belonging to the opening)
  # count parentheses
  pc <- 1
  cstart <- indexstart
  cend <- cstart
  while (pc != 0) {
    cend <- cend + 1
    if (substr(textnew,cend,cend) == '(') {
      pc <- pc+1
    } else {
      if (substr(textnew,cend,cend) == ')') {
        pc <- pc-1
      }
    }
  }
  indexend <- cend-1
  indexafter <- indexend+1
  # indexstart/indexend identify the content in the parentheses to be
  # processed and replaced
  textinside <- substr(textnew,indexstart,indexend)
  terms <- strexplodePC(textinside,separator=",",group="square")
  # We need now to make sure that the elements in the
  # vectors are separated using commata (otherwise big problem)!
  xtermstring <- strtrimM(terms[1])
  ytermstring <- strtrimM(terms[2])
  # remove parentheses
  xtermstring = substr(xtermstring,2,nchar(xtermstring)-1)
  ytermstring = substr(ytermstring,2,nchar(ytermstring)-1)
  # get single elements
  xtermelements = strexplodePC(xtermstring)
  ytermelements = strexplodePC(ytermstring)
  if (length(xtermelements) < 3) {
    stop('exchangeInterp0: the interp0 function requires at least 3 points on the x and y axis.')
  }
  if (length(xtermelements) != length(ytermelements)) {
    stop('exchangeInterp0: x and y arguments for interp0 function do not have same number of elements.');
  }
  # Construct piecewise statement
  pwText <- paste(ytermelements[1],",lt(",terms[3],",",xtermelements[1],"),",sep="")
  for (k in 2:length(xtermelements)-1) {
    pwText <- paste(pwText,"(",ytermelements[k],")",sep="")
    if (k<length(xtermelements)-1) {
      pwText <- paste(pwText,",andAZR(lt(",terms[3],",",xtermelements[k+1],"),ge(",terms[3],",",xtermelements[k],")),",sep="")
    }
  }
  pwText <- paste(pwText,",andAZR(lt(",terms[3],",",xtermelements[length(xtermelements)],"),ge(",terms[3],",",xtermelements[length(xtermelements)-1],")),(",ytermelements[length(ytermelements)],")",sep="")

  textnew <- paste(substr(text,1,indexstart-1), pwText, substr(text,indexend+1,nchar(text)))
  textnew <- strrepM(textnew,'interp0','piecewiseAZR')
  return(textnew)
}

##########################################################################
# HANDLE INTERP1 (LOOKUP TABLE W/ LINEAR INTERPOLATION)
##########################################################################
exchangeInterp1 <- function (text) {
  # Check if interp1 present
  if (is.null(strlocateall(text,"interp1(")$start)) {
    return(text)
  }
  # Check if present multiple times
  if (length(strlocateall(text,"interp1(")$start) > 1) {
    stop('exchangeInterp1: the "interp1" function is only allowed to be present once in each formula.')
  }

  # Present once => handle with piecewise
  textnew <- text

  # Get starting index of things inside
  indexstart <- strlocateall(text,"interp1(")$end+1

  # Cut out the content using the parentheses
  # Get the end index of the k-th statement (closing parentheses
  # belonging to the opening)
  # count parentheses
  pc <- 1
  cstart <- indexstart
  cend <- cstart
  while (pc != 0) {
    cend <- cend + 1
    if (substr(textnew,cend,cend) == '(') {
      pc <- pc+1
    } else {
      if (substr(textnew,cend,cend) == ')') {
        pc <- pc-1
      }
    }
  }
  indexend <- cend-1
  indexafter <- indexend+1
  # indexstart/indexend identify the content in the parentheses to be
  # processed and replaced
  textinside <- substr(textnew,indexstart,indexend)
  terms <- strexplodePC(textinside,separator=",",group="square")
  # We need now to make sure that the elements in the
  # vectors are separated using commata (otherwise big problem)!
  xtermstring <- strtrimM(terms[1])
  ytermstring <- strtrimM(terms[2])
  # remove parentheses
  xtermstring = substr(xtermstring,2,nchar(xtermstring)-1)
  ytermstring = substr(ytermstring,2,nchar(ytermstring)-1)
  # get single elements
  xtermelements = strexplodePC(xtermstring)
  ytermelements = strexplodePC(ytermstring)
  if (length(xtermelements) < 3) {
    stop('exchangeInterp1: the interp1 function requires at least 3 points on the x and y axis.')
  }
  if (length(xtermelements) != length(ytermelements)) {
    stop('exchangeInterp1: x and y arguments for interp1 function do not have same number of elements.');
  }
  # Construct piecewise statement
  pwText <- paste(ytermelements[1],",lt(",terms[3],",",xtermelements[1],"),",sep="")

  for (k in 2:length(xtermelements)-1) {
    pwText <- paste(pwText,"(",ytermelements[k+1],"-(",ytermelements[k],"))/(",xtermelements[k+1],"-(",xtermelements[k],"))*(",terms[3],"-(",xtermelements[k],"))+(",ytermelements[k],")",sep="")
    if (k<length(xtermelements)-1) {
      pwText <- paste(pwText,",andAZR(lt(",terms[3],",",xtermelements[k+1],"),ge(",terms[3],",",xtermelements[k],")),",sep="")
    }
  }
  pwText <- paste(pwText,",andAZR(lt(",terms[3],",",xtermelements[length(xtermelements)],"),ge(",terms[3],",",xtermelements[length(xtermelements)-1],")),(",ytermelements[length(ytermelements)],")",sep="")
  textnew <- paste(substr(text,1,indexstart-1), pwText, substr(text,indexend+1,nchar(text)))
  textnew <- strrepM(textnew,'interp1','piecewiseAZR')
  return(textnew)
}

##########################################################################
# HANDLE INTERPCS (LOOKUP TABLE W/ CUBIC SPLINE INTERPOLATION)
##########################################################################
exchangeInterpcs <- function (text) {
  # Check if interpcs present
  if (is.null(strlocateall(text,"interpcsAZR(")$start)) {
    return(text)
  }
  # Check if present multiple times
  if (length(strlocateall(text,"interpcsAZR(")$start) > 1) {
    stop('exchangeInterpcs: the "interpcs" function is only allowed to be present once in each formula.')
  }

  # Present once => handle with piecewise
  textnew <- text

  # Get starting index of things inside
  indexstart <- strlocateall(text,"interpcsAZR(")$end+1

  # Cut out the content using the parentheses
  # Get the end index of the k-th statement (closing parentheses
  # belonging to the opening)
  # count parentheses
  pc <- 1
  cstart <- indexstart
  cend <- cstart
  while (pc != 0) {
    cend <- cend + 1
    if (substr(textnew,cend,cend) == '(') {
      pc <- pc+1
    } else {
      if (substr(textnew,cend,cend) == ')') {
        pc <- pc-1
      }
    }
  }
  indexend <- cend-1
  indexafter <- indexend+1

  # indexstart/indexend identify the content in the parentheses to be
  # processed and replaced
  textinside <- substr(textnew,indexstart,indexend)
  terms <- strexplodePC(textinside,separator=",",group="square")
  # We need now to make sure that the elements in the
  # vectors are separated using commata (otherwise big problem)!
  xtermstring <- strtrimM(terms[1])
  ytermstring <- strtrimM(terms[2])
  # remove parentheses
  xtermstring = substr(xtermstring,2,nchar(xtermstring)-1)
  ytermstring = substr(ytermstring,2,nchar(ytermstring)-1)
  # get single elements
  xtermelements = strexplodePC(xtermstring)
  ytermelements = strexplodePC(ytermstring)
  if (length(xtermelements) < 3) {
    stop('exchangeInterpcs: the interpcs function requires at least 3 points on the x and y axis.')
  }
  if (length(xtermelements) != length(ytermelements)) {
    stop('exchangeInterpcs: x and y arguments for interpcs function do not have same number of elements.');
  }

  # Generate new expression for C-code interpcsAZR function
  newexpr <- paste(length(xtermelements),",",terms[3],sep="")
  newexpr <- paste(newexpr,',',xtermstring,',',ytermstring,sep="")
  textnew <- paste(substr(text,1,indexstart-1), newexpr, substr(text,indexend+1,nchar(text)),sep="")
  return(textnew)
}

##########################################################################
# HANDLE THE POWER OPERATOR
##########################################################################
# Some tests:
# ===========
# formula = "X^2"; convertPowerOperator(formula)
# formula = "X^2+"; convertPowerOperator(formula)
# formula = "X^2^2"; convertPowerOperator(formula)
# formula = "pow(X,2)^2"; convertPowerOperator(formula)
# formula = "exp(X,2)^2^(1-2)^1-2"; convertPowerOperator(formula)
# formula = "-EMAX*((Ac/Vc)^(hill-1))"; convertPowerOperator(formula)
# formula = "((Ac/Vc)^hill)"; convertPowerOperator(formula)
# formula = "(pow(EC50,hill)+(Ac/Vc)^hill)"; convertPowerOperator(formula)
# formula = "(EC50^hill+(Ac/Vc)^hill)"; convertPowerOperator(formula)
# formula = "-EMAX*((Ac/Vc)^(hill-1)*(hill*(1/Vc)))"; convertPowerOperator(formula)
# formula = "-EMAX*((Ac/Vc)^(hill-1)*(hill*(1/Vc)))/(EC50^hill+(Ac/Vc)^hill)"; convertPowerOperator(formula)
# formula = "-EMAX*pow((Ac/Vc),hill)*(pow((Ac/Vc),(hill-1))*(hill*(1/Vc)))/(pow(EC50,hill)+pow((Ac/Vc),hill))^2"; convertPowerOperator(formula)
# formula = "-EMAX*(Ac/Vc)^hill*((Ac/Vc)^(hill-1)*(hill*(1/Vc)))/(EC50^hill+(Ac/Vc)^hill)"; convertPowerOperator(formula)
# formula = "-EMAX*pow((Ac/Vc),hill)*(pow((Ac/Vc),(hill-1))*(hill*(1/Vc)))/(pow(EC50,hill)+pow((Ac/Vc),hill))^2"; convertPowerOperator(formula)
# formula = "(pow(EC50,hill)+(Ac/Vc)^hill)^2"; convertPowerOperator(formula)
# formula = "(EC50^hill+(Ac/Vc)^hill)^2"; convertPowerOperator(formula)
# formula = "-EMAX*(Ac/Vc)^hill*((Ac/Vc)^(hill-1)*(hill*(1/Vc)))/(EC50^hill+(Ac/Vc)^hill)^2"; convertPowerOperator(formula)
# formula <- "-(EMAX * ((Ac/Vc)^(hill - 1) * (hill * (1/Vc)))/(EC50^hill +     (Ac/Vc)^hill) - EMAX * (Ac/Vc)^hill * ((Ac/Vc)^(hill - 1) *     (hill * (1/Vc)))/(EC50^hill + (Ac/Vc)^hill)^2)"
# convertPowerOperator(formula)
convertPowerOperator <- function (formula) {
  # remove whitespaces from formula
  formula <- strremWhite(formula)

  # then do the more complicated stuff
  indices <- strlocateall(formula,'^')$start

  while (!is.null(indices)) {

    index <- indices[1]

    formula1 <- substr(formula,1,index-1)
    formula2 <- substr(formula,index+1,nchar(formula))

    # check formula1 from the right
    # search for first occurrence of +,-,*,/,( outside parentheses
    pc <- 0
    cend <- nchar(formula1)
    cstart <- cend
    run <- TRUE
    while (run) {
      if (pc==0) {
        if (substr(formula1,cstart,cstart) %in% c("+","-","*","/")) { cstart <- cstart+1; break }
        if (substr(formula1,cstart,cstart) == "(") {
          if (cstart == 1) { cstart <- cstart+1; break }
          if (substr(formula1,cstart-1,cstart-1) %in% c("+","-","*","/","(")) { cstart <- cstart+1; break }
        }
      }

      if (substr(formula1,cstart,cstart) == ')') pc <- pc+1
      if (substr(formula1,cstart,cstart) == '(') pc <- pc-1
      if (cstart<=1) break
      cstart <- cstart - 1
    }
    firstargument <- substr(formula1,cstart,cend)
    cendfirst <- cstart

    # check formula2 from the left
    # search for first occurrence of +,-,*,/,) outside parentheses
    pc <- 0
    cstart <- 1
    cend <- cstart
    run <- TRUE
    while (run) {
      if (pc==0 & substr(formula2,cend,cend) %in% c("+","-","*","/",")","^")) { cend <- cend-1; break }
      if (cend>=nchar(formula2)) break
      if (substr(formula2,cend,cend) == ')') pc <- pc+1
      if (substr(formula2,cend,cend) == '(') pc <- pc-1
      cend <- cend + 1
    }
    secondargument <- substr(formula2,cstart,cend)
    cstartsecond <- cend

    # construct power expression
    powerexp <- paste("pow(",firstargument,",",secondargument,")",sep="")

    # construct new formula
    formula <- paste(substr(formula1,1,cendfirst-1), powerexp, substr(formula2,cstartsecond+1,nchar(formula2)),sep="")

    # get new indices for '^' character
    indices <- strlocateall(formula,'^')$start
  }
  return(formula)
}
