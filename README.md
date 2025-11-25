# FuncCell: Function-Aware Single-Cell Embeddings

Create cell embeddings by combining single-cell RNA-seq expression with protein function embeddings from ProteinBERT.

**Project Members:** Andre Gala-Garza, Prema Immadisetty, Jason (Jingqiao) Zhao  
**Course:** BIOINF 593

## Overview

Raw gene expression counts lack functional context. FuncCell creates function-aware cell embeddings by weighting ProteinBERT protein embeddings by gene expression:

```
cell_embedding = Σ (normalized_expression[gene] × gene_embedding[gene])
```

This places cells in a "functional space" where similarity reflects shared protein functions, not just co-expression.

## Current Results

- **19,294** protein-coding genes embedded (512-dim via ProteinBERT)
- **99%** classification accuracy distinguishing healthy vs cancer breast cells (100 cells each, 5-fold CV)
- Clear separation in embedding space (distance between centroids: 0.86 in 512-dim)

## Project Structure

```
funcCell/
├── src/
│   ├── preprocess_data/
│   │   ├── config.py              # Configuration and paths
│   │   ├── census_query.py        # CELLxGENE Census queries
│   │   ├── gene_mapping.py        # Gene → protein sequence mapping
│   │   └── gene_list.py           # Gene list utilities
│   └── model/
│       ├── proteinbert_embeddings.py  # ProteinBERT embedding generation
│       └── cell_embeddings.py         # Cell-level aggregation
├── scripts/
│   ├── generate_embeddings.py     # Generate gene embeddings
│   └── create_cell_embeddings.py  # Generate cell embeddings
├── notebooks/
│   └── test_cell_embeddings.ipynb # Analysis notebook
├── data/                          # (gitignored)
│   ├── sequences/
│   │   ├── gene_to_sequence.pkl   # Gene → protein sequence dict
│   │   └── gene_list.txt          # List of mapped genes
│   ├── embeddings/
│   │   └── gene_to_embedding.pkl  # Gene → 512-dim embedding dict
│   └── analysis/
│       ├── cancer_direction_512d.npy   # Cancer direction vector
│       ├── healthy_centroid_512d.npy   # Healthy cell centroid
│       ├── cancer_centroid_512d.npy    # Cancer cell centroid
│       └── gene_contribution_analysis.csv  # Gene analysis results
└── models/
    └── proteinbert/               # ProteinBERT model weights
```

## Installation

```bash
# Clone and setup
git clone https://github.com/jason-zhao-git/FuncCell.git
cd FuncCell

# Install with uv
uv sync
source .venv/bin/activate
```

## Usage

### 1. Generate Gene Embeddings

```bash
python scripts/generate_embeddings.py
```

This:
- Loads protein sequences from `data/sequences/gene_to_sequence.pkl`
- Generates 512-dim ProteinBERT embeddings for each gene
- Saves to `data/embeddings/gene_to_embedding.pkl`

### 2. Create Cell Embeddings

```python
from src.model.cell_embeddings import CellEmbedder, load_gene_embeddings

# Load gene embeddings
gene_to_embedding = load_gene_embeddings("data/embeddings/gene_to_embedding.pkl")

# Create embedder
embedder = CellEmbedder(gene_to_embedding, embedding_dim=512)

# Generate cell embeddings from AnnData
cell_embeddings = embedder.create_cell_embeddings(adata)  # shape: (n_cells, 512)
```

### 3. Analyze (Interactive Notebook)

```bash
jupyter lab notebooks/test_cell_embeddings.ipynb
```

The notebook:
- Queries CELLxGENE Census for breast tissue data
- Generates cell embeddings
- Computes cancer direction vector in 512-dim space
- Identifies genes driving the cancer/healthy separation
- Exports gene lists for pathway analysis

## Method

### Cell Embedding Formula

For each cell with expression vector `e` and genes `G`:

```
weights[g] = e[g] / Σ e[g]  (normalize to sum to 1)
cell_emb = Σ weights[g] × gene_emb[g]
```

### Cancer Direction Analysis

```python
healthy_centroid = mean(healthy_cell_embeddings)  # (512,)
cancer_centroid = mean(cancer_cell_embeddings)    # (512,)
cancer_direction = cancer_centroid - healthy_centroid  # (512,)

# Score any cell
cancer_score = cell_embedding @ normalize(cancer_direction)
```

### Gene Contribution Score

For each gene:
- `alignment` = dot(gene_embedding, cancer_direction)
- `log2FC` = log2(cancer_expr / healthy_expr)
- `contribution_score` = alignment × log2FC

Genes with high positive scores drive cancer separation; high negative scores drive healthy.

## Configuration

Edit `src/preprocess_data/config.py`:

```python
CENSUS_QUERY_PARAMS = {
    "healthy_filter": "tissue_general == 'breast' and disease == 'normal' ...",
    "cancer_filter": "tissue_general == 'breast' and disease == 'breast cancer' ...",
    "min_genes": 200,   # QC: min genes per cell
    "min_cells": 300,   # QC: min cells per gene
}
```

## Dependencies

- `cellxgene-census` - Single-cell data access
- `scanpy`, `anndata` - Single-cell analysis
- `proteinbert`, `tensorflow` - Protein embeddings
- `numpy`, `pandas`, `scikit-learn` - Analysis

## License

Academic project for BIOINF 593.
