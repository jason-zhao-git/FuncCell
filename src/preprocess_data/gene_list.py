"""Load protein-coding genes from BioMart."""

import logging
import pandas as pd

from . import config

logger = logging.getLogger(__name__)


def get_protein_coding_genes() -> tuple[list, dict]:
    """
    Load all protein-coding gene symbols from BioMart export.

    Returns:
        Tuple of (gene_symbols, gene_descriptions)
        - gene_symbols: List of gene symbols for all protein-coding genes
        - gene_descriptions: Dict mapping gene symbol to gene description
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

    # Create mapping of gene symbols to descriptions
    gene_descriptions = {}
    for _, row in protein_coding.iterrows():
        gene_name = row['Gene name']
        gene_desc = row['Gene description']
        if pd.notna(gene_name) and pd.notna(gene_desc):
            gene_descriptions[gene_name] = gene_desc

    logger.info(f"Loaded {len(gene_symbols):,} protein-coding gene symbols from BioMart")
    logger.info(f"  - {len(gene_descriptions):,} genes have descriptions")

    return gene_symbols, gene_descriptions
