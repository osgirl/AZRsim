install.packages("DiagrammeR")
library(DiagrammeR)


get_replacement_vars <- function(model) {
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

grViz('
  digraph test {

      # several "node" statements
      node [shape = box,fontname = Helvetica]
      Ad; Cc;

      node [shape = circle,
      fixedsize = true,
      width = 0.9] // sets as circles
      input1; input2;

      // empty node to use to have arrows go to/from exclusively
      node [style=invis, width=0.1, height=0.1]
      empty;
      # several "edge" statements
      Ad -> Cc [label="Ka"]; input1 -> Ad;
      input2 -> Cc ; Cc -> empty [label="CL"];
      }
      ')

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
