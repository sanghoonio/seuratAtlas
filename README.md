## seuratAtlas

Interactive web-based visualizations for Seurat single-cell RNA-seq analysis results, using Apple's [Embedding Atlas](https://github.com/apple/embedding-atlas).

### Overview

`seuratAtlas` provides an interactive viewer for exploring Seurat objects through a web interface. The package launches a local web server that displays UMAP embeddings and metadata with real-time querying capabilities.

### Installation

You can install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("sanghoonio/seuratAtlas" subdir = "api")
```

### Quick Start

``` r
library(Seurat)
library(seuratAtlas)

# process your Seurat object with UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# launch the interactive viewer
seurat_atlas(pbmc)
```

This will:
1. Start a local web server at `http://localhost:3000`
2. Automatically open your default browser
3. Display an interactive interface for exploring your data

### Features

- **Interactive UMAP visualization**: Pan, zoom, and select cells
- **Real-time querying**: Fast SQL-based data filtering
- **Metadata exploration**: Visualize different metadata columns
- **WebSocket communication**: Responsive data updates

## Documentation

See the package vignette for a detailed tutorial:

``` r
vignette("pbmc-tutorial", package = "seuratAtlas")
```

Or view the function documentation:

``` r
?seurat_atlas
```

## Building the UI

The UI is prebuilt separately and copied into the package. To rebuild:

``` bash
cd ../ui
npm run build  # postbuild automatically copies dist folder to api/inst/www/
```
