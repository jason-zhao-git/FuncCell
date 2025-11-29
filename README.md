# FuncCell: Function-Aware Single-Cell Embeddings

Create cell embeddings by combining single-cell RNA-seq expression with protein function embeddings from ProteinBERT.

**Project Members:** Andre Gala-Garza, Prema Immadisetty, Jason (Jingqiao) Zhao
**Course:** BIOINF 593

## Overview

Raw gene expression counts lack functional context. FuncCell creates function-aware cell embeddings by weighting ProteinBERT protein embeddings by gene expression. We implement two approaches:

1. **Baseline (Weighted Sum):** Use expression as fixed attention weights
2. **Attention Pooling:** Learn attention weights from both gene identity and expression

This places cells in a "functional space" where similarity reflects shared protein functions, not just co-expression.

## Current Results

- **19,294** protein-coding genes embedded (512-dim via ProteinBERT)
- **20K training cells** (10K healthy + 10K cancer breast cells from CELLxGENE Census)
- **Attention Pooling** learns which genes matter regardless of expression level

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
│       ├── cell_embeddings.py         # Baseline cell-level aggregation
│       └── attention_pooling.py       # Attention pooling model
├── scripts/
│   ├── generate_embeddings.py     # Generate gene embeddings
│   └── create_cell_embeddings.py  # Generate cell embeddings
├── notebooks/
│   ├── test_cell_embeddings.ipynb   # Baseline analysis notebook
│   └── attention_pooling.ipynb      # Attention pooling training
├── data/                          # (gitignored)
│   ├── sequences/
│   │   ├── gene_to_sequence.pkl   # Gene → protein sequence dict
│   │   └── gene_list.txt          # List of mapped genes
│   ├── embeddings/
│   │   └── gene_to_embedding.pkl  # Gene → 512-dim embedding dict
│   └── analysis/
│       └── gene_contribution_analysis.csv  # Gene analysis results
├── models/
│   ├── proteinbert/               # ProteinBERT model weights
│   └── attention_pooling/         # Trained attention model weights
└── reference/                     # Technical documentation
    ├── embedding_generation.md
    ├── proteinbert_architecture.md
    ├── census_metadata_discovery.md
    └── gene_mapping_analysis.md
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

### 3. Train Attention Pooling

```bash
jupyter lab notebooks/attention_pooling.ipynb
```

The notebook:
- Queries CELLxGENE Census (or local SOMA) for breast tissue data
- Creates train/test datasets (10K cells each class)
- Trains AttentionPooling model with BCE loss
- Compares against baseline (weighted sum + LogReg)
- Visualizes embeddings with UMAP

### 4. Baseline Analysis

```bash
jupyter lab notebooks/test_cell_embeddings.ipynb
```

For simple weighted-sum analysis without learned attention.

## Method

### 1. Baseline: Weighted Sum

For each cell with expression vector `e` and genes `G`:

```
weights[g] = e[g] / Σ e[g]  (normalize to sum to 1)
cell_emb = Σ weights[g] × gene_emb[g]
```

Simple and fast, but assumes high expression = high importance.

### 2. Attention Pooling (Learned Weights)

Learn attention from both gene identity and expression:

```
attention[g] = f(gene_embedding[g], expression[g])
cell_emb = Σ softmax(attention)[g] × gene_emb[g]
```

**Architecture:**
```
Attention Network (shared across all ~19K genes):
  Input: [gene_embedding (512) | expression (1)] = 513 dims
  Linear(513 → 256) → ReLU → Linear(256 → 128) → ReLU → Linear(128 → 1)
  Output: attention score per gene

Classifier:
  Linear(512 → 64) → ReLU → Dropout(0.3) → Linear(64 → 1) → Sigmoid
```

**Parameters:** ~198K | **Loss:** Binary Cross-Entropy

**Why it helps:** Learns that some genes (e.g., tumor suppressors) matter even when lowly expressed, and housekeeping genes should be ignored even when highly expressed.

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

- `tiledbsoma`, `cellxgene-census` - Census/SOMA data access
- `scanpy`, `anndata` - Single-cell analysis
- `torch` - Attention pooling model
- `proteinbert`, `tensorflow` - Protein embeddings
- `numpy`, `pandas`, `scikit-learn` - Analysis

## License

Academic project for BIOINF 593.
