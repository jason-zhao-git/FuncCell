"""Load protein-coding genes from BioMart."""

import logging
import pandas as pd

from . import config

logger = logging.getLogger(__name__)


def get_protein_coding_genes() -> list:
    """
    Load all protein-coding gene symbols from BioMart export.

    Returns:
        List of gene symbols for all protein-coding genes
    """
    biomart_file = config.GENE_FILTER_PARAMS["biomart_file"]

    if not biomart_file.exists():
        raise FileNotFoundError(
            f"BioMart file not found: {biomart_file}\n"
            f"Please download protein-coding gene annotations from Ensembl BioMart."
        )

    logger.info(f"Loading protein-coding genes from BioMart...")

    # Read BioMart export (CSV format)
    df = pd.read_csv(biomart_file)

    # Filter to protein_coding genes
    protein_coding = df[df['Gene type'] == 'protein_coding']

    # Extract gene symbols (gene names)
    gene_symbols = protein_coding['Gene name'].dropna().unique().tolist()

    logger.info(f"Loaded {len(gene_symbols):,} protein-coding gene symbols from BioMart")

    return gene_symbols
