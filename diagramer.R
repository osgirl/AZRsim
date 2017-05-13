install.packages("DiagrammeR")
library(DiagrammeR)


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

model_states

model_chunks <- chunk_ode(model_states)

prepare_inputs <- function(model_states) {
  model_inputs <- purrr::set_names(
    stringr::str_extract_all(model_states, "input\\d+" ),
    names(model_states)
  )
  paste0(purrr::flatten_chr(
    purrr::map(names(model_inputs), ~ paste0(model_inputs[[.x]], " -> ", .x, ";"))
  ), collapse = " ")
}

get_cmpt_relationships <- function(model_states) {
  model_chunks <- chunk_ode(model_states)

  positive_relationships <-
    purrr::set_names(map(names(model_chunks), function(.x){
      chunks <- model_chunks[[.x]]
      positives <- which(names(chunks) == "p")
      chunks <- chunks[positives]
      to_remove <- which(str_detect(chunks, "input\\d+"))
      if (length(to_remove)) {
        chunks <- chunks[-to_remove]
      }
      return(chunks)
  }), names(model_chunks))
  negative_relationships <-
      purrr::set_names(map(names(model_chunks), function(.x){
      chunks <- model_chunks[[.x]]
      negatives <- which(names(chunks) == "n")
      chunks <- chunks[negatives]
      to_remove <- which(str_detect(chunks, "input\\d+"))
      if (length(to_remove)) {
        chunks <- chunks[-to_remove]
      }
      return(chunks)
  }), names(model_chunks))

  map(positive_relationships, function(.relationship){
    if (length(.relationship)) {
      map(.relationship, function(.chunk) {

      })
    }

  })
}

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
      }
      ", .model_states = sep_values(names(model_states)),
        .model_inputs = sep_values(get_inputs(model_states)),
        .inputs = prepare_inputs(model_states),
        .cmpt_relations = get_cmpt_relationships(model_states)
)

test_grviz
grViz(test_grviz)
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
