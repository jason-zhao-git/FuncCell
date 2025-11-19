"""Query CELLxGENE Census for breast cancer datasets."""

import logging
import pandas as pd
import cellxgene_census
import scanpy as sc
from anndata import AnnData
from pathlib import Path
from tqdm import tqdm

from . import config

logger = logging.getLogger(__name__)


def _load_protein_coding_genes() -> set:
    """
    Load protein-coding gene Ensembl IDs from BioMart export.

    Returns:
        Set of Ensembl gene IDs for protein-coding genes
    """
    biomart_file = config.GENE_FILTER_PARAMS["biomart_file"]

    if not biomart_file.exists():
        raise FileNotFoundError(
            f"BioMart file not found: {biomart_file}\n"
            f"Please download protein-coding gene annotations from Ensembl BioMart."
        )

    # Read BioMart export (CSV format)
    df = pd.read_csv(biomart_file)

    # Filter to protein_coding genes
    protein_coding = df[df['Gene type'] == 'protein_coding']

    # Extract Ensembl IDs
    ensembl_ids = set(protein_coding['Gene stable ID'].dropna())

    logger.info(f"Loaded {len(ensembl_ids):,} protein-coding genes from BioMart")

    return ensembl_ids


def query_census_data(
    dataset_type: str,
    output_path: Path,
) -> AnnData:
    """
    Query CELLxGENE Census for a specific dataset.

    Args:
        dataset_type: Either 'healthy' or 'tumor'
        output_path: Path to save the raw data

    Returns:
        AnnData object with raw counts
    """
    dataset_config = config.CENSUS_CONFIG[dataset_type]
    logger.info(f"Querying {dataset_config['name']}...")
    logger.info(f"Collection ID: {dataset_config['collection_id']}")

    # Open Census
    with cellxgene_census.open_soma(census_version="latest") as census:
        # Build organism filter
        organism = "Homo sapiens"

        # Build filters
        filters = dataset_config["filters"]
        filter_parts = []

        for key, value in filters.items():
            if isinstance(value, list):
                # Handle list values (e.g., disease in ['breast cancer', 'cancer'])
                # Use triple quotes to handle apostrophes in values safely
                value_str = ", ".join([f'"{v}"' for v in value])
                filter_parts.append(f'''{key} in [{value_str}]''')
            elif isinstance(value, bool):
                # Handle boolean values (e.g., is_primary_data == True)
                filter_parts.append(f'''{key} == {value}''')
            else:
                # Handle single string values - use triple quotes + double quotes
                # This safely handles apostrophes in values like "10x 3' v3"
                filter_parts.append(f'''{key} == "{value}"''')

        # Combine all filters
        obs_value_filter = " and ".join(filter_parts)
        logger.info(f"Filter: {obs_value_filter}")

        # Query Census for expression data (all genes - no var_value_filter)
        logger.info("Fetching expression data from Census...")
        adata = cellxgene_census.get_anndata(
            census=census,
            organism=organism,
            obs_value_filter=obs_value_filter,
        )

    logger.info(f"Initial data shape: {adata.shape} (cells × genes)")

    # Filter to protein-coding genes using BioMart annotations
    if config.GENE_FILTER_PARAMS["use_protein_coding_only"]:
        logger.info("Filtering to protein-coding genes using BioMart annotations...")
        protein_coding_genes = _load_protein_coding_genes()

        # Match by feature_id (Ensembl ID)
        gene_mask = adata.var.index.isin(protein_coding_genes)
        n_protein_coding = gene_mask.sum()
        logger.info(f"  - Found {n_protein_coding:,} / {len(protein_coding_genes):,} protein-coding genes in dataset")

        adata = adata[:, gene_mask].copy()
        logger.info(f"After protein-coding filter: {adata.shape} (cells × genes)")

    # Basic QC filtering
    logger.info("Applying basic QC filters...")
    sc.pp.filter_cells(adata, min_genes=config.QC_PARAMS["min_genes"])
    sc.pp.filter_genes(adata, min_cells=config.QC_PARAMS["min_cells"])

    logger.info(f"After QC: {adata.shape} (cells × genes)")
    logger.info(f"Dataset info:")
    logger.info(f"  - Total cells: {adata.n_obs:,}")
    logger.info(f"  - Total genes: {adata.n_vars:,}")
    logger.info(f"  - Median genes/cell: {adata.obs['n_genes'].median():.0f}")
    logger.info(f"  - Median UMI counts/cell: {adata.obs['total_counts'].median():.0f}")

    # Save raw data
    logger.info(f"Saving raw data to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)

    return adata


def query_all_datasets() -> tuple[AnnData, AnnData]:
    """
    Query both healthy and tumor datasets from CELLxGENE Census.

    Returns:
        Tuple of (healthy_adata, tumor_adata)
    """
    config.ensure_data_dirs()

    # Query healthy dataset
    logger.info("="*60)
    logger.info("Step 1/2: Querying healthy breast tissue data")
    logger.info("="*60)
    healthy_adata = query_census_data(
        dataset_type="healthy",
        output_path=config.OUTPUT_FILES["healthy_raw"]
    )

    # Query tumor dataset
    logger.info("\n" + "="*60)
    logger.info("Step 2/2: Querying breast cancer tumor data")
    logger.info("="*60)
    tumor_adata = query_census_data(
        dataset_type="tumor",
        output_path=config.OUTPUT_FILES["tumor_raw"]
    )

    logger.info("\n" + "="*60)
    logger.info("Census query completed successfully!")
    logger.info("="*60)
    logger.info(f"Healthy: {healthy_adata.n_obs:,} cells × {healthy_adata.n_vars:,} genes")
    logger.info(f"Tumor:   {tumor_adata.n_obs:,} cells × {tumor_adata.n_vars:,} genes")

    return healthy_adata, tumor_adata


if __name__ == "__main__":
    from utils.logging_utils import setup_logger

    setup_logger("census_query", level=logging.INFO)
    query_all_datasets()
