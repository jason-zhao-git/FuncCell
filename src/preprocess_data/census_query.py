"""Query CELLxGENE Census for breast tissue single-cell data."""

import logging
from pathlib import Path
from typing import Optional, Set

import cellxgene_census
import numpy as np
import scanpy as sc
from anndata import AnnData

from .config import CENSUS_QUERY_PARAMS, OUTPUT_FILES

logger = logging.getLogger(__name__)


class CensusQuerier:
    """Query and filter single-cell data from CELLxGENE Census."""

    def __init__(
        self,
        healthy_filter: Optional[str] = None,
        cancer_filter: Optional[str] = None,
        min_genes: int = 200,
        min_cells: int = 3,
    ):
        """
        Initialize Census querier.

        Args:
            healthy_filter: Value filter for healthy breast cells
            cancer_filter: Value filter for breast cancer cells
            min_genes: Minimum genes per cell (QC threshold)
            min_cells: Minimum cells per gene (QC threshold)
        """
        self.healthy_filter = healthy_filter or CENSUS_QUERY_PARAMS["healthy_filter"]
        self.cancer_filter = cancer_filter or CENSUS_QUERY_PARAMS["cancer_filter"]
        self.min_genes = min_genes
        self.min_cells = min_cells

    def query_census(
        self,
        value_filter: str,
        available_genes: Optional[Set[str]] = None,
    ) -> AnnData:
        """
        Query Census data with given filter.

        Args:
            value_filter: Census value filter string
            available_genes: Set of gene symbols with embeddings (for filtering)

        Returns:
            AnnData object with filtered cells and genes
        """
        logger.info(f"Querying Census with filter: {value_filter}")

        # Open Census and query
        with cellxgene_census.open_soma(census_version="stable") as census:
            adata = cellxgene_census.get_anndata(
                census=census,
                organism="Homo sapiens",
                measurement_name="RNA",
                X_name="raw",
                obs_value_filter=value_filter,
            )

        logger.info(f"Retrieved {adata.n_obs} cells × {adata.n_vars} genes")

        # Basic QC filtering
        logger.info("Applying QC filters...")
        sc.pp.filter_cells(adata, min_genes=self.min_genes)
        sc.pp.filter_genes(adata, min_cells=self.min_cells)
        logger.info(
            f"After QC: {adata.n_obs} cells × {adata.n_vars} genes "
            f"(min_genes={self.min_genes}, min_cells={self.min_cells})"
        )

        # Filter to genes with embeddings if provided
        if available_genes is not None:
            common_genes = [g for g in adata.var_names if g in available_genes]
            n_before = adata.n_vars
            adata = adata[:, common_genes]
            logger.info(
                f"Filtered to genes with embeddings: {n_before} → {adata.n_vars} genes "
                f"({len(common_genes) / n_before * 100:.1f}% coverage)"
            )

        return adata

    def fetch_healthy_data(
        self, available_genes: Optional[Set[str]] = None
    ) -> AnnData:
        """
        Fetch healthy breast tissue data.

        Args:
            available_genes: Set of gene symbols with embeddings

        Returns:
            AnnData with healthy breast cells
        """
        logger.info("=" * 80)
        logger.info("Fetching healthy breast tissue data")
        logger.info("=" * 80)

        adata = self.query_census(self.healthy_filter, available_genes)

        # Add metadata
        adata.obs["cancer_label"] = 0  # 0 = healthy
        adata.obs["tissue_type"] = "healthy_breast"

        logger.info(f"✓ Healthy data: {adata.n_obs} cells × {adata.n_vars} genes\n")
        return adata

    def fetch_cancer_data(self, available_genes: Optional[Set[str]] = None) -> AnnData:
        """
        Fetch breast cancer tissue data.

        Args:
            available_genes: Set of gene symbols with embeddings

        Returns:
            AnnData with breast cancer cells
        """
        logger.info("=" * 80)
        logger.info("Fetching breast cancer tissue data")
        logger.info("=" * 80)

        adata = self.query_census(self.cancer_filter, available_genes)

        # Add metadata
        adata.obs["cancer_label"] = 1  # 1 = cancer
        adata.obs["tissue_type"] = "breast_cancer"

        logger.info(f"✓ Cancer data: {adata.n_obs} cells × {adata.n_vars} genes\n")
        return adata

    def save_data(self, adata: AnnData, output_path: Path) -> None:
        """
        Save AnnData to h5ad file.

        Args:
            adata: AnnData object to save
            output_path: Path to save h5ad file
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(output_path)
        logger.info(f"Saved to {output_path}")
        logger.info(f"File size: {output_path.stat().st_size / 1024**2:.1f} MB")


def load_available_genes(gene_embedding_path: Path) -> Set[str]:
    """
    Load set of genes with available embeddings.

    Args:
        gene_embedding_path: Path to gene_to_embedding.pkl

    Returns:
        Set of gene symbols
    """
    import pickle

    logger.info(f"Loading available genes from {gene_embedding_path}")
    with open(gene_embedding_path, "rb") as f:
        gene_to_embedding = pickle.load(f)

    genes = set(gene_to_embedding.keys())
    logger.info(f"Loaded {len(genes)} genes with embeddings\n")
    return genes


def main():
    """Main function for testing Census queries."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    # Load available genes
    gene_embedding_path = OUTPUT_FILES["gene_to_embedding"]
    available_genes = load_available_genes(gene_embedding_path)

    # Initialize querier
    querier = CensusQuerier()

    # Fetch healthy data
    healthy_data = querier.fetch_healthy_data(available_genes)
    querier.save_data(healthy_data, OUTPUT_FILES["healthy_census"])

    # Fetch cancer data
    cancer_data = querier.fetch_cancer_data(available_genes)
    querier.save_data(cancer_data, OUTPUT_FILES["cancer_census"])

    logger.info("=" * 80)
    logger.info("SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Healthy cells: {healthy_data.n_obs:,}")
    logger.info(f"Cancer cells: {cancer_data.n_obs:,}")
    logger.info(f"Total cells: {healthy_data.n_obs + cancer_data.n_obs:,}")
    logger.info(f"Common genes: {healthy_data.n_vars}")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
