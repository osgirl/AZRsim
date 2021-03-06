#' determine the names of the inputs in the ode system of equations
#' @param model_states vector of ode model states
#' @examples
#' get_inputs(c("-ka*Ad + F11*input1", "ka*Ad - Cl*Ac/V + F22*input2"))
get_inputs <- function(model_states) {
  inputs <- unlist(stringr::str_extract_all(model_states, "input\\d+" ))
  return(unique(inputs))
}

#' light wrapper for collapse diagrammer separator requirements
#' @param .values vector to collapse to single string
#' @param .sep separator to value
sep_values <- function(.values, .sep = ";") {
  paste0(.values, sep = .sep, collapse = ' ')
}

#' get the position and assign it a sign value from a str_locate_all output
#' @param .x list of start and end positions for locate call
#' @param sign the sign value, (p/n)
#' @examples
#'  positives <- stringr::str_locate_all('-ka*Ad+F11*input1', '\\+')
#'  position_and_sign(positives)
position_and_sign <- function(.x, sign = "p") {
  if (nrow(.x)) {
    results <- c()
    # str_locate_all returns a matrix of start/end values, we
    # need to iterate by row, so could use apply, however this
    # is simple so will just loop
    for (r in 1:nrow(.x))  {
     results <- c(results, purrr::set_names(.x[r][[1]], sign))
    }
    return(results)
  }
  return(NULL)
}

#' determine the signs for each of the chunks in an ode
#' @param .string vector of strings defining ODEs
#' @examples
#' get_chunk_signs("-ka*Ad+F11*input1")
#' get_chunk_signs(list(
#'     Ad = "-ka*Ad+F11*input1",
#'     Ac = "ka*Ad + F22*input2 - CL*Ac/Vc")
#' )
get_chunk_signs <- function(.string) {
  positives <- stringr::str_locate_all(.string, "\\+")
  negatives <- stringr::str_locate_all(.string, "\\-")
  purrr::map2(positives, negatives, function(.x, .y) {
    positives <- position_and_sign(.x)
    negatives <- position_and_sign(.y, sign = "n")
    combined <- c(positives, negatives)
    combined <- combined[order(combined)]
    # if no 1st index it means that the start of the equation
    # has no +/-, and therefore is positive
    if (combined[[1]] != 1) {
      combined <- c(p = 1, combined)
    }
    return(combined)
  })
}

#' recombine equations in parens for one nested layer
#' @param .lv list of vectors of ode chunks
#' @examples
#' recombine_parens(list(c(p = "ka*Ad", p = "(1", n = "F11)*input1")))
recombine_parens <- function(.lv) {
  purrr::map(.lv, function(.v){
    to_combine <- c()
    # find if any start with ( and eventually combine it with
    # the next block to to turn, for example,
    # (1, F1)*input1 back to (1+F1)*input1
    for (i in seq_along(.v)) {
      if (substr(.v[i], 1, 1) == "(") {
        to_combine <- c(to_combine, i)
      }
    }
    if (length(to_combine)) {
      replacements <- purrr::map(to_combine, function(.c){
        indices <- .c:(.c+1)
        combination <- .v[indices]
        signs <- names(.v)[indices]
        sign_symbol <- ifelse(signs[2] == "p", "+", "-")
        # need to return list of result and sign as can't
        # set name inside a map function as won't bubble up
        # instead need to extract it after
        list(output = paste0(combination,
                             collapse = sign_symbol),
                     sign = signs[1])
      })
      replacements <- purrr::set_names(purrr::map_chr(replacements, ~ .x$output),
                               purrr::map_chr(replacements, ~ .x$sign))
      .v <- .v[-c(to_combine, to_combine+1)]
      .v <- c(.v, replacements)
    }
    return(.v)
  })
}

#' chunk an ode into parts split on +/- signs
#' @param .string vector of ode strings
#' @return list of chunked odes
chunk_ode <- function(.string) {
  splits <- stringr::str_split(.string, "[-|+]")
  signs <- get_chunk_signs(.string)
  results <- purrr::map2(splits, signs,
                         function(.splits, .signs){
    # if first split was negative the lhs will be empty
    # so want to remove the first empty block, eg
    # splitting -ka*Ad will have 2 blocks
    # a blank leading one, ka*Ad
    if (.splits[1] == "") {
      .splits <- .splits[-1]
    }
    if (!(length(.splits) == length(.signs))) {
      stop("signs and splits not matching")
    }
    purrr::set_names(stringr::str_trim(.splits), names(.signs))
  })
  results <- recombine_parens(results)
  if (!is.null(names(.string))) {
    return(purrr::set_names(results, names(.string)))
  }
  return(results)
}

#' prepare inputs from the system of equations for diagrammer
#' @param model_states named vector of ode system of equations
#' @details
#' nodes for diagrammer need to be specified, this will extract the
#' compartments and the inputs and generate the required syntax
prepare_inputs <- function(model_states) {
  model_inputs <- purrr::set_names(
    stringr::str_extract_all(model_states, "input\\d+" ),
    names(model_states)
  )
  paste0(purrr::flatten_chr(
    purrr::map(names(model_inputs), ~ paste0(model_inputs[[.x]], " -> ", .x, ";"))
  ), collapse = " ")
}
