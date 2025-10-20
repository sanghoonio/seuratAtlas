parse_params <- function(query_string) {
  if (is.null(query_string) || nchar(query_string) == 0) {
    return(list())
  }
  
  params <- list()
  pairs <- strsplit(query_string, '&')[[1]]
  for (pair in pairs) {
    kv <- strsplit(pair, '=')[[1]]
    if (length(kv) == 2) {
      params[[URLdecode(kv[1])]] <- utils::URLdecode(kv[2])
    }
  }
  params
}

read_post_body <- function(req) {
  if (req$REQUEST_METHOD == 'POST' && !is.null(req$rook.input)) {
    body <- rawToChar(req$rook.input$read())
    if (nchar(body) > 0) {
      return(jsonlite::fromJSON(body))
    }
  }
  return(NULL)
}
