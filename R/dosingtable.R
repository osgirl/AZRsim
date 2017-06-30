###############################################################################
###############################################################################
# Dosing table handling functions for the simulation
###############################################################################
###############################################################################

check_dosing_table <- function(dosing_table) {

  # Handle undefined dosing_table
  if (is.null(dosing_table)) return(NULL)

  # If dosing table contains ID this ID should be unique!
  if ("ID" %in% colnames(dosing_table)) {
    if (length(unique(dosing_table$ID))> 1)
      stop("check_dosing_table: dosing_table contains ID with non-unique entries. Plase check if you used the AZRsimulate instead of AZRsimpop function")
  }

  # Check if minimum required columns present
  # These are: TIME, INPUT, DOSE
  if (!("TIME" %in% names(dosing_table))) stop("TIME column required in a dosing table")
  if (!("INPUT" %in% names(dosing_table))) stop("INPUT column required in a dosing table")
  if (!("DOSE" %in% names(dosing_table))) stop("DOSE column required in a dosing table")

  # Sort the dosing table after TIME and INPUT
  dosing_table <- dplyr::arrange(dosing_table,TIME,INPUT)

  # Check if DURATION is present in dosing table
  # Handle DURATION and RATE columns
  if (!("DURATION" %in% names(dosing_table))) {
    # Check if RATE is present
    if ("RATE" %in% names(dosing_table)) {
      # RATE is present -> generate DURATION. 0 RATE indicates bolus (DURATION=0)
      # And remove RATE column afterwards
      DURATION <- dosing_table$DOSE/dosing_table$RATE
      DURATION[dosing_table$RATE==0] <- 0
      dosing_table$DURATION <- DURATION
      dosing_table$RATE <- NULL
    } else {
      # Neither RATE nor DURATION present => add DURATION with 0 entries
      dosing_table$DURATION <- 0
    }
  } else {
    # DURATION present
    # Check if RATE is present as well (not allowed)
    if ("RATE" %in% names(dosing_table)) stop("RATE and DURATION are not allowed to be present in a dosing table at the same time")
  }

  # Check if LAGTIME is present in dosing table
  if (!("LAGTIME" %in% names(dosing_table))) {
    # LAGTIME not present - add it with default "0" values
    dosing_table$LAGTIME <- 0
  }

  # Handle 0 DURATION TIME (set to 0.0001)
  dosing_table$DURATION[dosing_table$DURATION<.Machine$double.eps] <- 0.0001

  # Finally it needs to be checked that next dose comes after previous dosing is finished
  # Check for each input number independently ... allowing for example a very long infusion
  # with one input and daily dosing with the other.
  # (time+duration+lagtime < next dosing time)
  inputs <- unique(dosing_table$INPUT)
  for (k in 1:length(inputs)) {
    dosing_tableInputk <- dplyr::filter(dosing_table,INPUT==inputs[k])
    DoseTimes <- sort(unique(dosing_tableInputk$TIME))
    if (length(DoseTimes) > 1) {
      dosing_tableInputk$TEST_TIMING <- dosing_tableInputk$TIME+dosing_tableInputk$DURATION+dosing_tableInputk$LAGTIME
      for (k in 1:(length(DoseTimes)-1)) {
        if (max(dplyr::filter(dosing_tableInputk,TIME==DoseTimes[k])$TEST_TIMING) > DoseTimes[k+1])
          stop("check_dosing_table: dose administration of a dose happens after start of next dosing event")
      }
    }
  }

  # All adjusted and OK - return updated dosing_table
  return(dosing_table)
}
