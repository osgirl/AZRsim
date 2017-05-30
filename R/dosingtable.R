###############################################################################
###############################################################################
# Dosing table handling functions for the simulation
###############################################################################
###############################################################################

check_dosing_table <- function(dosingTable) {

  # Handle undefined dosingTable
  if (is.null(dosingTable)) return(NULL)

  # If dosing table contains ID this ID should be unique!
  if ("ID" %in% colnames(dosingTable)) {
    if (length(unique(dosingTable$ID))> 1)
      stop("check_dosing_table: dosingTable contains ID with non-unique entries. Plase check if you used the AZRsimulate instead of AZRsimpop function")
  }

  # Check if minimum required columns present
  # These are: TIME, INPUT, DOSE
  if (!("TIME" %in% names(dosingTable))) stop("TIME column required in a dosing table")
  if (!("INPUT" %in% names(dosingTable))) stop("INPUT column required in a dosing table")
  if (!("DOSE" %in% names(dosingTable))) stop("DOSE column required in a dosing table")

  # Sort the dosing table after TIME and INPUT
  dosingTable <- dplyr::arrange(dosingTable,TIME,INPUT)

  # Check if DURATION is present in dosing table
  # Handle DURATION and RATE columns
  if (!("DURATION" %in% names(dosingTable))) {
    # Check if RATE is present
    if ("RATE" %in% names(dosingTable)) {
      # RATE is present -> generate DURATION. 0 RATE indicates bolus (DURATION=0)
      # And remove RATE column afterwards
      DURATION <- dosingTable$DOSE/dosingTable$RATE
      DURATION[dosingTable$RATE==0] <- 0
      dosingTable$DURATION <- DURATION
      dosingTable$RATE <- NULL
    } else {
      # Neither RATE nor DURATION present => add DURATION with 0 entries
      dosingTable$DURATION <- 0
    }
  } else {
    # DURATION present
    # Check if RATE is present as well (not allowed)
    if ("RATE" %in% names(dosingTable)) stop("RATE and DURATION are not allowed to be present in a dosing table at the same time")
  }

  # Check if LAGTIME is present in dosing table
  if (!("LAGTIME" %in% names(dosingTable))) {
    # LAGTIME not present - add it with default "0" values
    dosingTable$LAGTIME <- 0
  }

  # Handle 0 DURATION TIME (set to 0.0001)
  dosingTable$DURATION[dosingTable$DURATION<.Machine$double.eps] <- 0.0001

  # Finally it needs to be checked that next dose comes after previous dosing is finished
  # Check for each input number independently ... allowing for example a very long infusion
  # with one input and daily dosing with the other.
  # (time+duration+lagtime < next dosing time)
  inputs <- unique(dosingTable$INPUT)
  for (k in 1:length(inputs)) {
    dosingTableInputk <- dplyr::filter(dosingTable,INPUT==inputs[k])
    DoseTimes <- sort(unique(dosingTableInputk$TIME))
    if (length(DoseTimes) > 1) {
      dosingTableInputk$TEST_TIMING <- dosingTableInputk$TIME+dosingTableInputk$DURATION+dosingTableInputk$LAGTIME
      for (k in 1:(length(DoseTimes)-1)) {
        if (max(dplyr::filter(dosingTableInputk,TIME==DoseTimes[k])$TEST_TIMING) > DoseTimes[k+1])
          stop("check_dosing_table: dose administration of a dose happens after start of next dosing event")
      }
    }
  }

  # All adjusted and OK - return updated dosingTable
  return(dosingTable)
}
