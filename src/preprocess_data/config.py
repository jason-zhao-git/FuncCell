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

# CELLxGENE Census dataset configuration
CENSUS_CONFIG = {
    "healthy": {
        "collection_id": "48259aa8-f168-4bf5-b797-af8e88da6637",
        "name": "Human breast cell atlas",
        "filters": {
            "tissue_general": "breast",     # Use tissue_general (not tissue)
            "disease": "normal",            # Healthy tissue
            "assay": "10x 3' v3",           # Standardize assay type for consistency
            "is_primary_data": True,        # Use only primary data, avoid duplicates
        },
    },
    "tumor": {
        "collection_id": "dea97145-f712-431c-a223-6b5f565f362a",
        "name": "Human breast cancer atlas",
        "filters": {
            "tissue_general": "breast",     # Use tissue_general (not tissue)
            "disease": "breast cancer",     # Generic breast cancer (641k cells)
            "assay": "10x 3' v3",           # Standardize assay type for consistency
            "is_primary_data": True,        # Use only primary data, avoid duplicates
        },
    },
}

# Quality control parameters
QC_PARAMS = {
    "min_genes": 200,  # Minimum genes per cell
    "min_cells": 3,    # Minimum cells per gene
    "max_cells": None,  # Maximum cells to query (None = all, set to e.g. 50000 for testing)
}

# Gene filtering parameters
GENE_FILTER_PARAMS = {
    # BioMart export with protein-coding gene annotations
    "biomart_file": DATA_DIR / "raw" / "mart_export.txt",
    "use_protein_coding_only": True,  # Filter to protein-coding genes after query
}

# HVG selection parameters
HVG_PARAMS = {
    "skip_hvg_selection": True,    # OPTION A: Skip HVG, use all protein-coding genes
    "n_top_genes": 2000,           # Number of highly variable genes (if not skipped)
    "flavor": "seurat_v3",         # HVG selection method
    "batch_key": None,             # No batch correction for now
}

# Normalization parameters
NORM_PARAMS = {
    "target_sum": 1e4,            # Target sum for normalization
}

# Gene mapping parameters
GENE_MAPPING_PARAMS = {
    "batch_size": 100,            # Genes per API request
    "max_retries": 3,             # Max retries for failed requests
    "timeout": 30,                # Request timeout in seconds
}

# Output file paths
OUTPUT_FILES = {
    "healthy_raw": RAW_DATA_DIR / "healthy_raw.h5ad",
    "tumor_raw": RAW_DATA_DIR / "tumor_raw.h5ad",
    "healthy_filtered": PROCESSED_DATA_DIR / "healthy_filtered.h5ad",
    "tumor_filtered": PROCESSED_DATA_DIR / "tumor_filtered.h5ad",
    "gene_list": SEQUENCES_DIR / "gene_list.txt",
    "gene_to_sequence": SEQUENCES_DIR / "gene_to_sequence.pkl",
    "metadata_summary": REPORTS_DIR / "metadata_summary.txt",
    "gene_mapping_report": REPORTS_DIR / "gene_mapping_report.txt",
}


def ensure_data_dirs():
    """Create all data directories if they don't exist."""
    for directory in [RAW_DATA_DIR, PROCESSED_DATA_DIR, SEQUENCES_DIR, REPORTS_DIR]:
        directory.mkdir(parents=True, exist_ok=True)
