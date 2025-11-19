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

# Output file paths
OUTPUT_FILES = {
    "gene_list": SEQUENCES_DIR / "gene_list.txt",
    "gene_to_sequence": SEQUENCES_DIR / "gene_to_sequence.pkl",
    "metadata_summary": REPORTS_DIR / "metadata_summary.txt",
    "gene_mapping_report": REPORTS_DIR / "gene_mapping_report.txt",
}


def ensure_data_dirs():
    """Create all data directories if they don't exist."""
    for directory in [RAW_DATA_DIR, PROCESSED_DATA_DIR, SEQUENCES_DIR, REPORTS_DIR]:
        directory.mkdir(parents=True, exist_ok=True)
