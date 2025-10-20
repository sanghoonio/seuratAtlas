#' Launch Embedding Atlas for Seurat Objects
#'
#' Opens an interactive web-based viewer for exploring Seurat single-cell
#' RNA-seq analysis results. The viewer displays UMAP embeddings and metadata
#' in a web interface accessible at http://localhost:3000/ui.
#'
#' @param seurat_obj A Seurat object with UMAP embeddings. The object must
#'   have a 'umap' reduction available (created with \code{RunUMAP()}).
#'
#' @return Starts a web server and opens the default browser. The function
#'   blocks until the server is stopped (Ctrl+C or Esc).
#'
#' @details
#' The function extracts UMAP coordinates and metadata from the Seurat object,
#' stores them in an in-memory DuckDB database, and launches a web server with
#' WebSocket support for real-time querying.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(seuratAtlas)
#'
#' # Load and process your data
#' pbmc <- RunUMAP(pbmc, dims = 1:10)
#'
#' # Launch the interactive viewer
#' seurat_atlas(pbmc)
#' }
#'
#' @export
seurat_atlas <- function(seurat_obj) {
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ':memory:')
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  
  coords <- as.data.frame(Seurat::Embeddings(seurat_obj, reduction = 'umap'))
  metadata <- seurat_obj@meta.data
  
  data <- data.frame(
    cell_id = rownames(coords),
    projection_x = coords[, 1],
    projection_y = coords[, 2],
    metadata,
    check.names = FALSE
  )
  
  DBI::dbWriteTable(con, 'seurat_obj', data, overwrite = TRUE)
  
  cache <- cachem::cache_mem()
  create_server(con, cache, run_ui = TRUE)
}
