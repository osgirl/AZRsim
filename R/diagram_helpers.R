#' get the position and assign it a sign value from a str_locate_all output
#' @param .x list of start and end positions for locate call
#' @param sign the sign value, (p/n)
#' @examples
#'  positives <- str_locate_all('-ka*Ad+F11*input1', '\\+')
#'  position_and_sign(positives)
position_and_sign <- function(.x, sign = "p") {
  if (nrow(.x)) {
    results <- c()
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
  positives <- str_locate_all(.string, "\\+")
  negatives <- str_locate_all(.string, "\\-")
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
  map(.lv, function(.v){
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
      replacements <- map(to_combine, function(.c){
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
      replacements <- purrr::set_names(map_chr(replacements, ~ .x$output),
                               map_chr(replacements, ~ .x$sign))
      .v <- .v[-c(to_combine, to_combine+1)]
      .v <- c(.v, replacements)
    }
    return(.v)
  })
}

#' chunk an ode into parts split on +/- signs
#' @param .string vector of ode strings
chunk_ode <- function(.string) {
  splits <- str_split(.string, "[-|+]")
  signs <- get_chunk_signs(.string)
  purrr::map2(splits, signs, function(.splits, .signs){
    # if first split was negative the lhs will be empty
    if (.splits[1] == "") {
      .splits <- .splits[-1]
    }
    if (!(length(.splits) == length(.signs))) {
      stop("signs and splits not matching")
    }
    purrr::set_names(.splits, names(.signs))
  })
}
