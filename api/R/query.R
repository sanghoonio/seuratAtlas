get_key <- function(sql, command) {
  paste0(digest::digest(sql, algo = 'sha256'), '_', command)
}

retrieve <- function(cache, query, get) {
  key <- get_key(query$sql, query$type)
  result <- cache$get(key)

  if (cachem::is.key_missing(result)) { # not in cache
    result <- get(query$sql)
    if (isTRUE(query$persist)) {
      cache$set(key, result)
    }
  }
  result
}

get_json <- function(con, sql) {
  result <- DBI::dbGetQuery(con, sql)
  result[is.na(result)] <- ''
  as.character(jsonlite::toJSON(result, auto_unbox = TRUE))
}

get_arrow_bytes <- function(con, sql) {
  result <- DBI::dbGetQuery(con, sql)
  result[is.na(result)] <- ''
  arrow::write_to_raw(result, format = 'stream')
}

handle_query <- function(con, cache, query) {
  sql <- query$sql
  command <- query$type
  
  result <- tryCatch({
    if (command == 'exec') {
      DBI::dbExecute(con, sql)
      list(type = 'done')
    } else if (command == 'json') {
      json <- retrieve(cache, query, function(s) get_json(con, s))
      list(type = 'json', data = json)
    } else if (command == 'arrow') {
      arrow_bytes <- retrieve(cache, query, function(s) get_arrow_bytes(con, s))
      list(type = 'arrow', data = arrow_bytes)
    } else {
      stop(paste('Unknown command: ', command))
    }
  }, error = function(e) {
    list(type = 'error', data = e$message)
  })
  result
}
