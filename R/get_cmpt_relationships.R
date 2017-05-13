
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
  generate_mappings <- function(.r) {
    map(.r, function(.relationship){
      if (length(.relationship) && !is.na(.relationship)) {
        has_relationship <- map_chr(.relationship, function(.chunk) {
          if (is.na(.chunk)) {
            return(NA_character_)
          }
          relationship_index <- which(str_detect(.chunk, names(model_chunks)))
          names(model_chunks)[relationship_index]
        })
        return(has_relationship[!is.na(has_relationship)])
      }
    return(NA_character_)
  })
  }
  # negative relationships should be taken out as already accounted
  # for via positive flow
  negative_relationships <- map(negative_relationships, function(.chunks) {
    map_chr(.chunks, function(.chunk) {
      if (.chunk %in% unlist(positive_relationships)) {
        return(NA_character_)
      }
      return(.chunk)
    })
  })
  positive_mappings <- generate_mappings(positive_relationships)
  negative_mappings <- generate_mappings(negative_relationships)
  positive_mappings <- positive_mappings[!is.na(positive_mappings)]
  negative_mappings <- negative_mappings[!is.na(negative_mappings)]
  dplyr::bind_rows(
    map2_df(positive_mappings, names(positive_mappings), function(.x, .y) {
      map_df(.x, function(.x) {
        tibble::data_frame(from = .x, to = .y)
      })
  }),
    map_df(negative_mappings, function(.x) {
      # negative mappings are only going to exist if values
      # are exiting from its own compartment, (hopefully),
      # therefore should be piped to empty
      tibble::data_frame(from = .x, to = "empty")
  })
  )
}
