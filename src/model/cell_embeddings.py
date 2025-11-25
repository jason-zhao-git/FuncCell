"""Generate cell-level embeddings from gene embeddings and expression data."""

import logging
from pathlib import Path
from typing import Dict

import numpy as np
from anndata import AnnData
from scipy.sparse import issparse

logger = logging.getLogger(__name__)


class CellEmbedder:
    """Aggregate gene embeddings to create cell-level representations."""

    def __init__(self, gene_to_embedding: Dict[str, np.ndarray], embedding_dim: int = 512):
        """
        Initialize cell embedder.

        Args:
            gene_to_embedding: Dictionary mapping gene symbols to embeddings
            embedding_dim: Dimension of gene embeddings (default: 512)
        """
        self.gene_to_embedding = gene_to_embedding
        self.embedding_dim = embedding_dim
        logger.info(f"Initialized CellEmbedder with {len(gene_to_embedding)} genes")

    def create_cell_embeddings(
        self,
        adata: AnnData,
        batch_size: int = 1000,
    ) -> np.ndarray:
        """
        Create cell-level embeddings via expression-weighted sum of gene embeddings.

        Formula: cell_embedding = normalized_expression @ gene_embeddings
        where normalized_expression is expression / sum(expression) per cell

        Args:
            adata: AnnData with expression data (cells × genes)
            batch_size: Number of cells to process per batch

        Returns:
            Array of shape (n_cells, embedding_dim) with cell embeddings
        """
        # Filter to genes with embeddings
        common_genes = [g for g in adata.var_names if g in self.gene_to_embedding]
        n_common = len(common_genes)
        logger.info(
            f"Found {n_common}/{adata.n_vars} genes with embeddings "
            f"({n_common / adata.n_vars * 100:.1f}% coverage)"
        )

        if n_common == 0:
            raise ValueError("No common genes found between data and embeddings")

        # Filter AnnData to common genes
        adata_filtered = adata[:, common_genes]

        # Get expression matrix: (n_cells, n_genes)
        expression = adata_filtered.X
        if issparse(expression):
            logger.info("Converting sparse matrix to dense for aggregation")
            expression = expression.toarray()

        # Stack gene embeddings: (n_genes, embedding_dim)
        gene_embeddings = np.array([self.gene_to_embedding[g] for g in common_genes])
        logger.info(f"Gene embeddings shape: {gene_embeddings.shape}")
        logger.info(f"Expression matrix shape: {expression.shape}")

        # Process in batches for memory efficiency
        n_cells = expression.shape[0]
        n_batches = (n_cells + batch_size - 1) // batch_size
        cell_embeddings = np.zeros((n_cells, self.embedding_dim), dtype=np.float32)

        logger.info(f"Processing {n_cells} cells in {n_batches} batches...")

        for batch_idx in range(n_batches):
            start_idx = batch_idx * batch_size
            end_idx = min(start_idx + batch_size, n_cells)
            batch_expression = expression[start_idx:end_idx]

            # Normalize expression per cell: sum to 1
            row_sums = batch_expression.sum(axis=1, keepdims=True)
            # Add small epsilon to avoid division by zero
            normalized_expression = batch_expression / (row_sums + 1e-10)

            # Weighted sum: (batch_size, n_genes) @ (n_genes, embedding_dim)
            #             = (batch_size, embedding_dim)
            batch_embeddings = normalized_expression @ gene_embeddings

            cell_embeddings[start_idx:end_idx] = batch_embeddings

            if (batch_idx + 1) % 10 == 0 or batch_idx == n_batches - 1:
                logger.info(
                    f"  Batch {batch_idx + 1}/{n_batches}: "
                    f"processed {end_idx}/{n_cells} cells"
                )

        logger.info(f"✓ Generated cell embeddings: {cell_embeddings.shape}")
        return cell_embeddings

    def validate_embeddings(self, embeddings: np.ndarray) -> Dict[str, float]:
        """
        Validate cell embeddings and compute statistics.

        Args:
            embeddings: Cell embeddings array (n_cells, embedding_dim)

        Returns:
            Dictionary with validation statistics
        """
        stats = {}

        # Check for invalid values
        stats["n_cells"] = embeddings.shape[0]
        stats["embedding_dim"] = embeddings.shape[1]
        stats["has_nan"] = np.isnan(embeddings).any()
        stats["has_inf"] = np.isinf(embeddings).any()
        stats["n_nan"] = np.isnan(embeddings).sum()
        stats["n_inf"] = np.isinf(embeddings).sum()

        # Compute magnitude statistics
        norms = np.linalg.norm(embeddings, axis=1)
        stats["norm_mean"] = float(norms.mean())
        stats["norm_std"] = float(norms.std())
        stats["norm_min"] = float(norms.min())
        stats["norm_max"] = float(norms.max())

        # Log validation results
        logger.info("=" * 80)
        logger.info("VALIDATION RESULTS")
        logger.info("=" * 80)
        logger.info(f"  Shape: {embeddings.shape}")
        logger.info(f"  NaN values: {stats['n_nan']}")
        logger.info(f"  Inf values: {stats['n_inf']}")
        logger.info(f"  Embedding magnitude (mean ± std): {stats['norm_mean']:.3f} ± {stats['norm_std']:.3f}")
        logger.info(f"  Embedding magnitude range: [{stats['norm_min']:.3f}, {stats['norm_max']:.3f}]")

        if stats["has_nan"] or stats["has_inf"]:
            logger.warning("⚠ Embeddings contain NaN or Inf values!")
        else:
            logger.info("✓ All embeddings are valid (no NaN/Inf)")

        logger.info("=" * 80)

        return stats


def load_gene_embeddings(gene_embedding_path: Path) -> Dict[str, np.ndarray]:
    """
    Load gene embeddings from pickle file.

    Args:
        gene_embedding_path: Path to gene_to_embedding.pkl

    Returns:
        Dictionary mapping gene symbols to embeddings
    """
    import pickle

    logger.info(f"Loading gene embeddings from {gene_embedding_path}")
    with open(gene_embedding_path, "rb") as f:
        gene_to_embedding = pickle.load(f)

    logger.info(f"Loaded {len(gene_to_embedding)} gene embeddings")

    # Validate embedding dimension
    first_embedding = next(iter(gene_to_embedding.values()))
    embedding_dim = first_embedding.shape[0]
    logger.info(f"Embedding dimension: {embedding_dim}")

    return gene_to_embedding


def save_cell_embeddings(embeddings: np.ndarray, output_path: Path) -> None:
    """
    Save cell embeddings to numpy file.

    Args:
        embeddings: Cell embeddings array
        output_path: Path to save .npy file
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(output_path, embeddings)
    logger.info(f"Saved cell embeddings to {output_path}")
    logger.info(f"File size: {output_path.stat().st_size / 1024**2:.1f} MB")
