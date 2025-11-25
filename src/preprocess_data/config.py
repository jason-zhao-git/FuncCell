"""Configuration for data preprocessing pipeline."""

from pathlib import Path

# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent

# Data directories
DATA_DIR = PROJECT_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
SEQUENCES_DIR = DATA_DIR / "sequences"
REPORTS_DIR = DATA_DIR / "reports"
CENSUS_DIR = DATA_DIR / "census"
CELL_EMBEDDINGS_DIR = DATA_DIR / "cell_embeddings"
EMBEDDINGS_DIR = DATA_DIR / "embeddings"

# Gene filtering parameters
GENE_FILTER_PARAMS = {
    # BioMart export with protein-coding gene annotations
    "biomart_file": DATA_DIR / "raw" / "mart_export.txt",
}

# Gene mapping parameters
GENE_MAPPING_PARAMS = {
    "batch_size": 100,            # Genes per API request
    "max_retries": 3,             # Max retries for failed requests
    "timeout": 30,                # Request timeout in seconds
}

# Census query parameters
CENSUS_QUERY_PARAMS = {
    # Filter for healthy breast tissue
    "healthy_filter": "tissue_general == 'breast' and disease == 'normal' and assay == \"10x 3' v3\" and is_primary_data == True",

    # Filter for breast cancer tissue
    "cancer_filter": "tissue_general == 'breast' and disease == 'breast cancer' and assay == \"10x 3' v3\" and is_primary_data == True",

    # Quality control thresholds
    "min_genes": 200,             # Minimum genes per cell
    "min_cells": 300,             # Minimum cells per gene

    # Processing parameters
    "batch_size": 1000,           # Cells to process per batch
}

# Output file paths
OUTPUT_FILES = {
    # Gene sequences and mappings
    "gene_list": SEQUENCES_DIR / "gene_list.txt",
    "gene_to_sequence": SEQUENCES_DIR / "gene_to_sequence.pkl",
    "gene_to_embedding": EMBEDDINGS_DIR / "gene_to_embedding.pkl",

    # Reports
    "metadata_summary": REPORTS_DIR / "metadata_summary.txt",
    "gene_mapping_report": REPORTS_DIR / "gene_mapping_report.txt",
    "embedding_generation_report": EMBEDDINGS_DIR / "embedding_generation_report.txt",

    # Census data
    "healthy_census": CENSUS_DIR / "healthy_breast.h5ad",
    "cancer_census": CENSUS_DIR / "cancer_breast.h5ad",

    # Cell embeddings
    "healthy_cell_embeddings": CELL_EMBEDDINGS_DIR / "healthy_cells.npy",
    "cancer_cell_embeddings": CELL_EMBEDDINGS_DIR / "cancer_cells.npy",
    "cell_embeddings_report": CELL_EMBEDDINGS_DIR / "cell_embeddings_report.txt",
}


def ensure_data_dirs():
    """Create all data directories if they don't exist."""
    for directory in [
        RAW_DATA_DIR,
        PROCESSED_DATA_DIR,
        SEQUENCES_DIR,
        REPORTS_DIR,
        CENSUS_DIR,
        CELL_EMBEDDINGS_DIR,
        EMBEDDINGS_DIR,
    ]:
        directory.mkdir(parents=True, exist_ok=True)
