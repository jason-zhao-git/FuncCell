# funcCell: Single-Cell Cancer Analysis with Protein Function Embeddings

A bioinformatics pipeline for creating function-aware embeddings of single cells by combining transcriptomic profiles with protein functional information from ProteinBERT.

**Project Members:** Andre Gala-Garza, Prema Immadisetty, Jason (Jingqiao) Zhao
**Course:** BIOINF 593
**Presentation Date:** December 02, 2025

## Project Overview

This project addresses the challenge that raw gene expression counts lack direct functional context. We create function-aware cell embeddings by fusing single-cell RNA-seq data with ProteinBERT's pre-trained protein function representations.

**Goal:** Distinguish cancerous from healthy breast cells in a learned "functional space" where cancer and healthy cells form distinct, separable clusters.

### Datasets

- **Tumor Data:** [A single-cell and spatially resolved atlas of human breast cancers](https://cellxgene.cziscience.com/collections/dea97145-f712-431c-a223-6b5f565f362a)
- **Healthy Data:** [Human breast cell atlas](https://cellxgene.cziscience.com/collections/48259aa8-f168-4bf5-b797-af8e88da6637)
- **Protein Model:** [ProteinBERT](https://huggingface.co/GrimSqueaker/proteinBERT)

## Project Structure

```
funcCell/
├── src/
│   ├── preprocess_data/        # Data querying & gene mapping
│   │   ├── config.py           # Configuration parameters
│   │   ├── census_query.py     # CELLxGENE Census data fetching
│   │   ├── hvg_selection.py    # HVG selection on combined data
│   │   └── gene_mapping.py     # Gene → protein sequence mapping
│   └── model/                  # Embedding models
│       └── cell_embeddings.py  # Cell-level embedding aggregation
├── scripts/
│   ├── run_data_pipeline.py    # Data preprocessing pipeline
│   └── create_cell_embeddings.py  # Cell embedding generation
├── notebooks/
│   └── test_cell_embeddings.ipynb  # Interactive analysis notebook
├── data/                       # Generated data (gitignored)
│   ├── raw/                    # Raw queried data
│   ├── processed/              # Filtered data
│   ├── sequences/              # Gene-to-protein mappings
│   ├── embeddings/             # Gene embeddings (512-dim)
│   ├── cell_embeddings/        # Cell embeddings
│   └── analysis/               # Cancer direction analysis
├── reference/                  # Project documentation
└── pyproject.toml              # Project dependencies
```

## Installation

### Prerequisites

- Python ≥ 3.9
- [uv](https://github.com/astral-sh/uv) package manager

### Setup Environment

```bash
# Install dependencies using uv
uv sync

# Activate the virtual environment
source .venv/bin/activate  # On Linux/Mac
# or
.venv\Scripts\activate     # On Windows
```

### Manual Installation (if uv is not available)

```bash
pip install -r requirements.txt
```

## Usage

### Day 1: Data Preprocessing Pipeline

Run the complete data preprocessing pipeline:

```bash
python scripts/run_data_pipeline.py
```

This script will:

1. **Query CELLxGENE Census** for healthy and tumor breast datasets
   - Filter to **10x 3' v3 assay** for technical consistency
   - Filter to **primary data only** to avoid duplicates
   - Filter to **protein-coding genes** (~15,000-20,000 genes)
   - Apply QC filters (min_genes=200, min_cells=3)
   - Save raw data to `data/raw/`

2. **Select Highly Variable Genes (HVGs)**
   - Combine healthy + tumor datasets
   - Normalize and log-transform
   - Calculate top 2,000 HVGs on combined data
   - Save filtered data to `data/processed/`

3. **Map Genes to Protein Sequences**
   - Use MyGene.info API to map gene symbols → UniProt IDs
   - Fetch protein sequences from UniProt
   - Expected coverage: >85% of HVGs
   - Save mappings to `data/sequences/`

**Expected Runtime:** 30-60 minutes (depends on network speed)

### Output Files

After running the pipeline, you'll find:

```
data/
├── raw/
│   ├── healthy_raw.h5ad              # Raw healthy cells (post-QC)
│   └── tumor_raw.h5ad                # Raw tumor cells (post-QC)
├── processed/
│   ├── healthy_filtered.h5ad         # Healthy cells × HVGs
│   └── tumor_filtered.h5ad           # Tumor cells × HVGs
├── sequences/
│   ├── gene_to_sequence.pkl          # {gene: protein_sequence} dict
│   └── gene_list.txt                 # Final list of mapped genes
└── reports/
    ├── metadata_summary.txt          # Overall pipeline statistics
    ├── gene_mapping_report.txt       # Gene mapping details
    └── pipeline_YYYYMMDD_HHMMSS.log  # Full execution log
```

### Running Individual Steps

You can also run individual pipeline steps:

```python
# Step 1: Query Census data
from src.preprocess_data.census_query import query_all_datasets
healthy_adata, tumor_adata = query_all_datasets()

# Step 2: Select HVGs
from src.preprocess_data.hvg_selection import select_hvgs
healthy_filtered, tumor_filtered, hvg_list = select_hvgs(healthy_adata, tumor_adata)

# Step 3: Map genes to proteins
from src.preprocess_data.gene_mapping import map_genes_to_proteins
gene_to_sequence, stats = map_genes_to_proteins(hvg_list)
```

## Key Features

### Data Pipeline (Day 1)

- ✅ Automated querying from CELLxGENE Census
- ✅ Protein-coding gene filtering using BioMart annotations (21,490 genes)
- ✅ **Option A Active:** Skip HVG selection, use all protein-coding genes (~15k-18k)
- ✅ Robust gene-to-protein mapping with retry logic
- ✅ Comprehensive logging and error handling
- ✅ Detailed statistics and QC reports

### Gene Embeddings (Day 2)

- ✅ ProteinBERT embedding generation (512-dim per gene)
- ✅ 19,294 protein-coding genes embedded

### Cell Embeddings (Day 3)

- ✅ Expression-weighted aggregation of gene embeddings
- ✅ Cell embedding = Σ (normalized_expression × gene_embedding)
- ✅ CELLxGENE Census integration for breast tissue data
- ✅ Cancer direction vector analysis in 512-dim space
- ✅ 99% classification accuracy (Logistic Regression, 5-fold CV)

### Remaining Tasks

- [ ] Triplet loss training for improved separation
- [ ] Extended evaluation across multiple cancer types
- [ ] Visualization dashboard

## Technical Approach

### Data Quality & Gene Selection Strategy

**Quality Filters (applied at query time):**
1. **Assay standardization**: Filter to `10x 3' v3` only
   - Reduces technical batch effects
   - Ensures comparable gene detection across samples
   - Prevents model from learning assay-specific artifacts
2. **Primary data only**: Filter to `is_primary_data == True`
   - Avoids duplicate/reprocessed data
   - Ensures original count matrices

**Gene Selection Pipeline:**
1. Query Census without `var_value_filter` (Census lacks `feature_biotype`)
2. Filter to protein-coding genes using BioMart (`data/raw/mart_export.txt`)
3. **Option A (Active):** Use all ~15k-18k protein-coding genes detected
4. **Option B (Disabled):** Calculate top 2,000 HVGs on combined data
   - Set `skip_hvg_selection: False` in config to enable
5. Final output: ~13k-16k genes with protein sequences

### Cell Embedding Method

1. **Gene Embeddings:** ProteinBERT generates 512-dim embeddings for each gene's protein sequence
2. **Expression Weighting:** Normalize expression per cell (weights sum to 1)
3. **Aggregation:** Cell embedding = Σ (normalized_expression[gene] × gene_embedding[gene])
4. **Analysis:** Cancer direction vector computed as (cancer_centroid - healthy_centroid) in 512-dim space

## Configuration

Edit `src/preprocess_data/config.py` to modify:

- Dataset collection IDs
- QC parameters (min_genes, min_cells)
- HVG selection parameters (n_top_genes, method)
- Gene mapping parameters (batch_size, retries)
- Output paths

## Dependencies

Key packages:

- `cellxgene-census` - CELLxGENE data access
- `scanpy` - Single-cell analysis
- `anndata` - Annotated data structures
- `mygene` - Gene annotation API
- `biopython` - Sequence handling
- `pandas`, `numpy` - Data manipulation
- `tqdm` - Progress bars

## Troubleshooting

### Census Query Issues

If Census queries fail or are slow:

```python
# Check Census status
import cellxgene_census
cellxgene_census.get_census_version_description("latest")

# Try specific version instead of "latest"
```

### Gene Mapping Issues

If gene mapping has low success rate (<85%):

- Check internet connectivity
- Verify gene symbols are standard HGNC symbols
- Check `data/reports/gene_mapping_report.txt` for missing genes

### Memory Issues

If you encounter memory errors:

- Reduce the number of cells by adding stricter QC filters
- Process datasets separately instead of combined HVG selection

## Contributing

This is a course project. For questions, contact the project members.

## License

Academic project for BIOINF 593.

## Acknowledgments

- CELLxGENE Census for providing standardized single-cell data
- ProteinBERT authors for the protein function model
- Course instructors and teaching assistants
