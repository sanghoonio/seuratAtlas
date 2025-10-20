create_server <- function(con, cache, run_ui) {
  app <- list(
    onWSOpen = function(ws) { # websocket server
      ws$onMessage(function(binary, message) {
        tryCatch({
          query <- jsonlite::fromJSON(message)
          result <- handle_query(con, cache, query)
          
          if (result$type == 'done') {
            ws$send(as.character(jsonlite::toJSON(list(), auto_unbox = TRUE)))
          } else if (result$type == 'json') {
            ws$send(result$data)
          } else if (result$type == 'arrow') {
            ws$send(result$data)
          } else if (result$type == 'error') {
            ws$send(as.character(jsonlite::toJSON(list(error = result$data), auto_unbox = TRUE)))
          }
        }, error = function(e) {
          print(e)
          ws$send(as.character(jsonlite::toJSON(list(error = paste('Server error:', e$message)), auto_unbox = TRUE)))
        })
      })
      
      ws$onClose(function() {
      })
    },
    call = function(req) { # http server
      headers <- list(
        'Access-Control-Allow-Origin' = '*',
        'Access-Control-Request-Method' = '*',
        'Access-Control-Allow-Methods' = 'OPTIONS, POST, GET',
        'Access-Control-Allow-Headers' = '*',
        'Access-Control-Max-Age' = '2592000'
      )
      
      if (req$REQUEST_METHOD == 'OPTIONS') {
        return(list(status = 200L, headers = headers, body = ''))
      }
      
      path <- req$PATH_INFO
      method <- req$REQUEST_METHOD
      params <- parse_params(req$QUERY_STRING)
      
      if (path == '/api/echo') {
        if (req$REQUEST_METHOD == 'GET') {
          query_string <- req$QUERY_STRING
          msg <- NULL
          
          if (!is.null(query_string) && grepl('msg=', query_string)) {
            msg_part <- sub('.*msg=([^&]*).*', '\\1', query_string)
            msg <- URLdecode(msg_part)
          }
          
          if (is.null(msg) || msg == '') {
            return(list(status = 400L, headers = headers, 
                        body = jsonlite::toJSON(list(error = 'Missing msg parameter'))))
          }
          
          headers$`Content-Type` <- 'application/json'
          return(list(
            status = 200L, 
            headers = headers, 
            body = jsonlite::toJSON(list(message = msg), auto_unbox = TRUE)
          ))
        } else {
          list(
            status = 405L,
            headers = list('Content-Type' = 'text/plain'),
            body = 'Method Not Allowed'
          )
        }
      }
      else if (path == '/ui/' || grepl('/ui/', path)) { # serve ui index
        www_dir <- system.file('www', package = 'seuratAtlas')
        index_path <- if (nchar(www_dir) > 0) {
          file.path(www_dir, 'index.html')
        } else {
          file.path('..', 'ui', 'dist', 'index.html')  # fallback for development
        }

        if (file.exists(index_path)) {
          index_content <- readBin(index_path, 'raw', n = file.info(index_path)$size)
          list(
            status = 200L,
            headers = list('Content-Type' = 'text/html'),
            body = index_content
          )
        } else {
          list(status = 404L, body = 'UI assets not found. Follow the README in the /ui directory to build the frontend.')
        }
      }
    },
    staticPaths = list(
      '/ui/assets' = {
        www_dir <- system.file('www', package = 'seuratAtlas')
        if (nchar(www_dir) > 0) {
          file.path(www_dir, 'assets')
        } else {
          file.path('..', 'ui', 'dist', 'assets')  # fallback for development
        }
      }
    )
  )
  
  # cat('Server listening at ws://localhost:3000 and http://localhost:3000\n')
  cat('Navigating to http://localhost:3000/ui/...\n')
  if (isTRUE(run_ui)) utils::browseURL('http://localhost:3000/ui/')
  httpuv::runServer('0.0.0.0', 3000, app)
}
