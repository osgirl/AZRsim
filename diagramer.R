library(DiagrammeR)
library(AZRsim)
model <- create_model("modelPKtest.txt")
model <- create_model("modelPK_3cpt_linear.txt")
model_states <- AZRsim:::get_all_states(model)$stateODEs

get_states_replacements <- function(model) {
  model_states <- AZRsim:::get_all_states(model)$stateODEs
  replacements <- AZRsim:::get_all_variables(model)$varformulas
  # remove input and outputs
  removals <- stringr::str_detect(
    tolower(names(replacements)),
    pattern = "[input|output]"
    )
  replacements <- replacements[-which(removals)]
  if (!length(replacements)) {
    return(model_states)
  }
  inverted <- purrr::set_names(names(replacements), replacements)
  stringr::str_replace_all(model_states, inverted)
}


get_inputs <- function(model_states) {
  inputs <- unlist(stringr::str_extract_all(model_states, "input\\d+" ))
  return(unique(inputs))
}

sep_values <- function(.values, .sep = ";") {
  paste0(.values, sep = .sep, collapse = ' ')
}


model_chunks <- chunk_ode(model_states)

prepare_inputs <- function(model_states) {
  model_inputs <- unlist(purrr::set_names(
    stringr::str_extract_all(model_states, "input\\d+" ),
    names(model_states)
  ) )
  paste0(purrr::flatten_chr(
    purrr::map(names(model_inputs), ~ paste0(model_inputs[[.x]], " -> ", .x, ";"))
  ), collapse = " ")
}

cmpt_relations <- get_cmpt_relationships(model_states) %>%
  mutate(result = paste0(from, "->", to, ";"))
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
grViz(test_grviz)

# hand coded example
grViz("
digraph test {

      # several 'node' statements
      node [shape = box,fontname = Helvetica]
      Ad; Cc;

      node [shape = circle,
      fixedsize = true,
      width = 0.9] // sets as circles
      input1; input2;

      // empty node to use to have arrows go to/from exclusively
      node [style=invis, width=0.1, height=0.1]
      empty;
      # several 'edge' statements
      Ad -> Cc [label='Ka']; input1 -> Ad;
      input2 -> Cc ; Cc -> empty [label='CL'];
      }
      ")
