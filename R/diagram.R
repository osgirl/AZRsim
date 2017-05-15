#' generate the digraph syntax from a model
#' @param model azrmodel
#' @export
generate_diagram_syntax <- function(model) {
  model_states <- AZRsim:::get_all_states(model)$stateODEs
  cmpt_relations <- get_cmpt_relationships(model_states) %>%
    dplyr::mutate(result = paste0(from, "->", to, ";"))
    paste0(cmpt_relations$result, collapse = " ")

  test_grviz <- infuser::infuse("
    digraph test {

        # several 'node' statements
        node [shape = box,fontname = Helvetica]
        {{.model_states}}

        node [shape = circle,
        fixedsize = true,
        width = 0.9]
        {{.model_inputs}}

        node [style=invis, width=0.1, height=0.1]
        empty;

        # inputs
        {{.inputs}}
        # several 'edge' statements
        {{.cmpt_relations}}
        }
        ", .model_states = sep_values(names(model_states)),
          .model_inputs = sep_values(get_inputs(model_states)),
          .inputs = prepare_inputs(model_states),
          .cmpt_relations =  paste0(cmpt_relations$result, collapse = " ")
  )

  test_grviz
}

#' create a diagram of the ODE system of a model
#' @param model azrmodel
#' @export
diagram <- function(model) {
  diagram <- generate_diagram_syntax(model)
  DiagrammeR::grViz(diagram)
}
