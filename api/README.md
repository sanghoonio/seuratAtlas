## API

This is the R package directory.

### Dependencies

- `Seurat`
- `httpuv`
- `duckdb`
- `arrow`
- `jsonlite`
- `digest`
- `cachem`
- `DBI`

### Architecture

This R package runs an httpuv server that:
- Handles WebSocket connections for real-time data transfer
- Uses an in-memory DuckDB database for fast querying
- Serves a React-based frontend from `inst/www/`
