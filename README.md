## Seurat Atlas

Interactive visualizations for Seurat single-cell analyses, through Apple's [Embedding Atlas](https://github.com/apple/embedding-atlas).

### Overview

`seuratAtlas` provides an interactive web interface for exploring Seurat objects. The package launches a local web server that displays Seurat embeddings and metadata with real-time querying capabilities and interactions.

### Installation

Install the development version from GitHub:

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

# launch the embedding atlas
seurat_atlas(pbmc)
```

This will:
1. Start a local web server at `http://localhost:3000`
2. Automatically open your default browser
3. Display an interactive interface for exploring your data

<img width="3024" height="1890" alt="atlas" src="https://github.com/user-attachments/assets/37e82859-d8b8-4310-a972-d4eb664bc52f" />

### Features

- **Interactive UMAP visualization**: Pan, zoom, and select cells
- **Real-time querying**: Fast SQL-based data filtering with websocket
- **Metadata exploration**: Query metadata columns with dynamic selections

## Documentation

View the package vignette for a detailed tutorial:

``` r
vignette("pbmc-tutorial", package = "seuratAtlas")
```

## Building the UI

The UI is prebuilt separately and copied into the package. To rebuild:

``` bash
cd ../ui
npm run build  # postbuild automatically copies dist folder to api/inst/www/
```
