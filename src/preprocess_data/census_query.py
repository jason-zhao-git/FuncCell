"""Query CELLxGENE Census for breast cancer datasets."""

import logging
import cellxgene_census
import scanpy as sc
from anndata import AnnData
from pathlib import Path
from tqdm import tqdm

from . import config

logger = logging.getLogger(__name__)


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
                value_str = ", ".join([f"'{v}'" for v in value])
                filter_parts.append(f"{key} in [{value_str}]")
            elif isinstance(value, bool):
                # Handle boolean values (e.g., is_primary_data == True)
                filter_parts.append(f"{key} == {value}")
            else:
                # Handle single string values
                filter_parts.append(f"{key} == '{value}'")

        # Combine all filters
        obs_value_filter = " and ".join(filter_parts)
        logger.info(f"Filter: {obs_value_filter}")

        # Query Census for expression data
        logger.info("Fetching expression data from Census...")
        adata = cellxgene_census.get_anndata(
            census=census,
            organism=organism,
            obs_value_filter=obs_value_filter,
            var_value_filter=f"feature_biotype == '{config.GENE_FILTER_PARAMS['feature_biotype']}'",
        )

    logger.info(f"Initial data shape: {adata.shape} (cells × genes)")

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
