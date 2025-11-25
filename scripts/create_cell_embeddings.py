#!/usr/bin/env python3
"""Generate cell-level embeddings from gene embeddings and Census data.

This script:
1. Loads gene embeddings (19,294 genes × 512-dim)
2. Queries CELLxGENE Census for breast tissue data (healthy + cancer)
3. Aggregates gene embeddings to cell level via expression-weighted sum
4. Saves cell embeddings and generates validation report
"""

import logging
import sys
from datetime import datetime
from pathlib import Path
from time import time

import numpy as np

# Add src to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from model.cell_embeddings import (
    CellEmbedder,
    load_gene_embeddings,
    save_cell_embeddings,
)
from preprocess_data.census_query import CensusQuerier, load_available_genes
from preprocess_data.config import OUTPUT_FILES, CENSUS_QUERY_PARAMS, ensure_data_dirs

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def print_header(title: str) -> None:
    """Print a formatted header."""
    logger.info("\n" + "=" * 80)
    logger.info(f"{title:^80}")
    logger.info("=" * 80 + "\n")


def generate_validation_report(
    healthy_embeddings: np.ndarray,
    cancer_embeddings: np.ndarray,
    healthy_stats: dict,
    cancer_stats: dict,
    execution_time: float,
    output_path: Path,
) -> None:
    """Generate and save a validation report."""
    report_lines = [
        "Cell Embedding Generation Report",
        "=" * 80,
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Execution time: {execution_time / 60:.1f} minutes",
        "",
        "INPUT DATA",
        "-" * 80,
        "",
        f"  Gene embeddings:              {OUTPUT_FILES['gene_to_embedding']}",
        f"  Healthy Census data:          {OUTPUT_FILES['healthy_census']}",
        f"  Cancer Census data:           {OUTPUT_FILES['cancer_census']}",
        "",
        "CONFIGURATION",
        "-" * 80,
        "",
        f"  Embedding dimension:          {healthy_embeddings.shape[1]}",
        f"  Batch size:                   {CENSUS_QUERY_PARAMS['batch_size']}",
        f"  Min genes per cell:           {CENSUS_QUERY_PARAMS['min_genes']}",
        f"  Min cells per gene:           {CENSUS_QUERY_PARAMS['min_cells']}",
        "",
        "RESULTS - HEALTHY CELLS",
        "-" * 80,
        "",
        f"  Total cells:                  {healthy_stats['n_cells']:,}",
        f"  Embedding dimension:          {healthy_stats['embedding_dim']}",
        f"  NaN values:                   {healthy_stats['n_nan']}",
        f"  Inf values:                   {healthy_stats['n_inf']}",
        f"  Magnitude (mean ± std):       {healthy_stats['norm_mean']:.3f} ± {healthy_stats['norm_std']:.3f}",
        f"  Magnitude range:              [{healthy_stats['norm_min']:.3f}, {healthy_stats['norm_max']:.3f}]",
        "",
        "RESULTS - CANCER CELLS",
        "-" * 80,
        "",
        f"  Total cells:                  {cancer_stats['n_cells']:,}",
        f"  Embedding dimension:          {cancer_stats['embedding_dim']}",
        f"  NaN values:                   {cancer_stats['n_nan']}",
        f"  Inf values:                   {cancer_stats['n_inf']}",
        f"  Magnitude (mean ± std):       {cancer_stats['norm_mean']:.3f} ± {cancer_stats['norm_std']:.3f}",
        f"  Magnitude range:              [{cancer_stats['norm_min']:.3f}, {cancer_stats['norm_max']:.3f}]",
        "",
        "SUMMARY",
        "-" * 80,
        "",
        f"  Total cells processed:        {healthy_stats['n_cells'] + cancer_stats['n_cells']:,}",
        f"  Healthy:                      {healthy_stats['n_cells']:,}",
        f"  Cancer:                       {cancer_stats['n_cells']:,}",
        f"  All embeddings valid:         {'✓' if not (healthy_stats['has_nan'] or healthy_stats['has_inf'] or cancer_stats['has_nan'] or cancer_stats['has_inf']) else '✗'}",
        "",
        "OUTPUT FILES",
        "-" * 80,
        "",
        f"  Healthy embeddings:           {OUTPUT_FILES['healthy_cell_embeddings']}",
        f"    - Shape:                    {healthy_embeddings.shape}",
        f"    - Size:                     {OUTPUT_FILES['healthy_cell_embeddings'].stat().st_size / 1024**2:.1f} MB",
        "",
        f"  Cancer embeddings:            {OUTPUT_FILES['cancer_cell_embeddings']}",
        f"    - Shape:                    {cancer_embeddings.shape}",
        f"    - Size:                     {OUTPUT_FILES['cancer_cell_embeddings'].stat().st_size / 1024**2:.1f} MB",
        "",
        "NEXT STEPS",
        "-" * 80,
        "",
        "1. Load embeddings:",
        "   healthy_emb = np.load('data/cell_embeddings/healthy_cells.npy')",
        "   cancer_emb = np.load('data/cell_embeddings/cancer_cells.npy')",
        "",
        "2. Create labels:",
        "   y_healthy = np.zeros(len(healthy_emb))",
        "   y_cancer = np.ones(len(cancer_emb))",
        "",
        "3. Train classifier:",
        "   X = np.vstack([healthy_emb, cancer_emb])",
        "   y = np.hstack([y_healthy, y_cancer])",
        "   # Train logistic regression, random forest, etc.",
        "",
        "=" * 80,
    ]

    # Write report
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("\n".join(report_lines))

    logger.info(f"\n✓ Validation report saved to {output_path}")


def main():
    """Main execution function."""
    start_time = time()

    print_header("Cell-Level Embedding Generation")

    # Ensure directories exist
    ensure_data_dirs()

    # Step 1: Load gene embeddings
    print_header("Step 1: Loading Gene Embeddings")
    gene_to_embedding = load_gene_embeddings(OUTPUT_FILES["gene_to_embedding"])
    available_genes = set(gene_to_embedding.keys())

    # Step 2: Query Census data
    print_header("Step 2: Querying CELLxGENE Census")
    querier = CensusQuerier()

    # Fetch or load healthy data
    if OUTPUT_FILES["healthy_census"].exists():
        logger.info(f"Loading cached healthy data from {OUTPUT_FILES['healthy_census']}")
        import scanpy as sc
        healthy_data = sc.read_h5ad(OUTPUT_FILES["healthy_census"])
    else:
        healthy_data = querier.fetch_healthy_data(available_genes)
        querier.save_data(healthy_data, OUTPUT_FILES["healthy_census"])

    # Fetch or load cancer data
    if OUTPUT_FILES["cancer_census"].exists():
        logger.info(f"Loading cached cancer data from {OUTPUT_FILES['cancer_census']}")
        import scanpy as sc
        cancer_data = sc.read_h5ad(OUTPUT_FILES["cancer_census"])
    else:
        cancer_data = querier.fetch_cancer_data(available_genes)
        querier.save_data(cancer_data, OUTPUT_FILES["cancer_census"])

    # Step 3: Create cell embeddings
    print_header("Step 3: Creating Cell Embeddings")
    embedder = CellEmbedder(gene_to_embedding, embedding_dim=512)

    logger.info("Processing healthy cells...")
    healthy_embeddings = embedder.create_cell_embeddings(
        healthy_data, batch_size=CENSUS_QUERY_PARAMS["batch_size"]
    )
    save_cell_embeddings(healthy_embeddings, OUTPUT_FILES["healthy_cell_embeddings"])

    logger.info("\nProcessing cancer cells...")
    cancer_embeddings = embedder.create_cell_embeddings(
        cancer_data, batch_size=CENSUS_QUERY_PARAMS["batch_size"]
    )
    save_cell_embeddings(cancer_embeddings, OUTPUT_FILES["cancer_cell_embeddings"])

    # Step 4: Validate embeddings
    print_header("Step 4: Validating Embeddings")

    logger.info("\nHealthy cells validation:")
    healthy_stats = embedder.validate_embeddings(healthy_embeddings)

    logger.info("\nCancer cells validation:")
    cancer_stats = embedder.validate_embeddings(cancer_embeddings)

    # Step 5: Generate report
    print_header("Step 5: Generating Validation Report")
    execution_time = time() - start_time

    generate_validation_report(
        healthy_embeddings,
        cancer_embeddings,
        healthy_stats,
        cancer_stats,
        execution_time,
        OUTPUT_FILES["cell_embeddings_report"],
    )

    # Final summary
    print_header("COMPLETE")
    logger.info(f"Total execution time: {execution_time / 60:.1f} minutes")
    logger.info(f"Healthy cells: {healthy_stats['n_cells']:,}")
    logger.info(f"Cancer cells: {cancer_stats['n_cells']:,}")
    logger.info(f"Total cells: {healthy_stats['n_cells'] + cancer_stats['n_cells']:,}")
    logger.info(f"Embedding dimension: {healthy_embeddings.shape[1]}")
    logger.info("\n✓ Cell embeddings ready for cancer classification training!")


if __name__ == "__main__":
    main()
