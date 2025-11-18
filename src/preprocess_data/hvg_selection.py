"""Highly Variable Gene (HVG) selection on combined datasets."""

import logging
import scanpy as sc
from anndata import AnnData
import numpy as np
from pathlib import Path

from . import config

logger = logging.getLogger(__name__)


def select_hvgs(
    healthy_adata: AnnData,
    tumor_adata: AnnData,
) -> tuple[AnnData, AnnData, list]:
    """
    Select highly variable genes from combined healthy + tumor data.

    This ensures we capture cancer-specific signals by computing HVGs
    on the combined dataset, rather than on each dataset separately.

    Args:
        healthy_adata: Healthy cells AnnData object
        tumor_adata: Tumor cells AnnData object

    Returns:
        Tuple of (filtered_healthy_adata, filtered_tumor_adata, hvg_list)
    """
    logger.info("="*60)
    logger.info("Selecting Highly Variable Genes (HVGs)")
    logger.info("="*60)

    # Add labels to track cell origin
    healthy_adata.obs['cell_type'] = 'healthy'
    tumor_adata.obs['cell_type'] = 'tumor'

    # Concatenate datasets
    logger.info("Concatenating healthy and tumor datasets...")
    logger.info(f"  Healthy: {healthy_adata.n_obs:,} cells")
    logger.info(f"  Tumor:   {tumor_adata.n_obs:,} cells")

    combined_adata = healthy_adata.concatenate(
        tumor_adata,
        batch_key='dataset',
        batch_categories=['healthy', 'tumor'],
        index_unique=None,
    )

    logger.info(f"  Combined: {combined_adata.n_obs:,} cells × {combined_adata.n_vars:,} genes")

    # Normalize and log-transform
    logger.info("Normalizing combined dataset...")
    sc.pp.normalize_total(combined_adata, target_sum=config.NORM_PARAMS["target_sum"])
    sc.pp.log1p(combined_adata)

    # Calculate highly variable genes on combined data
    logger.info(f"Calculating top {config.HVG_PARAMS['n_top_genes']} HVGs...")
    sc.pp.highly_variable_genes(
        combined_adata,
        n_top_genes=config.HVG_PARAMS["n_top_genes"],
        flavor=config.HVG_PARAMS["flavor"],
        batch_key=config.HVG_PARAMS["batch_key"],
    )

    # Get HVG list
    hvg_genes = combined_adata.var_names[combined_adata.var['highly_variable']].tolist()
    n_hvgs = len(hvg_genes)

    logger.info(f"Selected {n_hvgs} highly variable genes")
    logger.info(f"  Mean variance ratio: {combined_adata.var.loc[hvg_genes, 'variances_norm'].mean():.3f}")

    # Filter original datasets to HVGs
    logger.info("\nFiltering datasets to HVG set...")

    # Normalize and log-transform original datasets
    sc.pp.normalize_total(healthy_adata, target_sum=config.NORM_PARAMS["target_sum"])
    sc.pp.log1p(healthy_adata)
    sc.pp.normalize_total(tumor_adata, target_sum=config.NORM_PARAMS["target_sum"])
    sc.pp.log1p(tumor_adata)

    # Filter to HVGs
    healthy_filtered = healthy_adata[:, hvg_genes].copy()
    tumor_filtered = tumor_adata[:, hvg_genes].copy()

    logger.info(f"  Healthy filtered: {healthy_filtered.n_obs:,} cells × {healthy_filtered.n_vars:,} genes")
    logger.info(f"  Tumor filtered:   {tumor_filtered.n_obs:,} cells × {tumor_filtered.n_vars:,} genes")

    # Save filtered datasets
    logger.info("\nSaving filtered datasets...")
    healthy_filtered.write_h5ad(config.OUTPUT_FILES["healthy_filtered"])
    logger.info(f"  Saved: {config.OUTPUT_FILES['healthy_filtered']}")

    tumor_filtered.write_h5ad(config.OUTPUT_FILES["tumor_filtered"])
    logger.info(f"  Saved: {config.OUTPUT_FILES['tumor_filtered']}")

    # Save gene list
    with open(config.OUTPUT_FILES["gene_list"], 'w') as f:
        for gene in hvg_genes:
            f.write(f"{gene}\n")
    logger.info(f"  Saved gene list: {config.OUTPUT_FILES['gene_list']}")

    logger.info("\nHVG selection completed successfully!")

    return healthy_filtered, tumor_filtered, hvg_genes


def load_filtered_datasets() -> tuple[AnnData, AnnData, list]:
    """
    Load previously filtered datasets and gene list.

    Returns:
        Tuple of (healthy_adata, tumor_adata, gene_list)
    """
    logger.info("Loading filtered datasets...")

    healthy_adata = sc.read_h5ad(config.OUTPUT_FILES["healthy_filtered"])
    tumor_adata = sc.read_h5ad(config.OUTPUT_FILES["tumor_filtered"])

    with open(config.OUTPUT_FILES["gene_list"], 'r') as f:
        gene_list = [line.strip() for line in f if line.strip()]

    logger.info(f"  Healthy: {healthy_adata.n_obs:,} cells × {healthy_adata.n_vars:,} genes")
    logger.info(f"  Tumor:   {tumor_adata.n_obs:,} cells × {tumor_adata.n_vars:,} genes")
    logger.info(f"  Genes:   {len(gene_list):,}")

    return healthy_adata, tumor_adata, gene_list


if __name__ == "__main__":
    from utils.logging_utils import setup_logger
    import scanpy as sc

    setup_logger("hvg_selection", level=logging.INFO)

    # Load raw data
    logger.info("Loading raw datasets...")
    healthy_adata = sc.read_h5ad(config.OUTPUT_FILES["healthy_raw"])
    tumor_adata = sc.read_h5ad(config.OUTPUT_FILES["tumor_raw"])

    # Select HVGs
    select_hvgs(healthy_adata, tumor_adata)
