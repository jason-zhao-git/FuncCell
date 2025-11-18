#!/usr/bin/env python3
"""
Main script to run the complete Day 1 data preprocessing pipeline.

This script orchestrates:
1. Querying CELLxGENE Census for healthy and tumor breast datasets
2. Selecting highly variable genes (HVGs) from combined data
3. Mapping genes to protein sequences via MyGene.info and UniProt
4. Generating comprehensive reports

Usage:
    python scripts/run_data_pipeline.py
"""

import sys
import logging
from pathlib import Path
from datetime import datetime

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from src.preprocess_data import config
from src.preprocess_data.census_query import query_all_datasets
from src.preprocess_data.hvg_selection import select_hvgs
from src.preprocess_data.gene_mapping import map_genes_to_proteins
from utils.logging_utils import setup_logger, log_step
from utils.file_utils import save_text

# Setup logging
logger = setup_logger(
    "data_pipeline",
    log_file=config.REPORTS_DIR / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log",
    level=logging.INFO
)


@log_step(logger, "Complete Data Pipeline")
def run_pipeline():
    """Run the complete data preprocessing pipeline."""

    # Ensure output directories exist
    config.ensure_data_dirs()

    logger.info("\n")
    logger.info("╔" + "═" * 78 + "╗")
    logger.info("║" + " " * 20 + "funcCell Day 1 Data Pipeline" + " " * 30 + "║")
    logger.info("║" + " " * 15 + "Single-Cell Cancer Analysis with Protein Embeddings" + " " * 12 + "║")
    logger.info("╚" + "═" * 78 + "╝")
    logger.info("\n")

    # Track overall statistics
    stats = {}
    start_time = datetime.now()

    # ========================================================================
    # Step 1: Query CELLxGENE Census
    # ========================================================================
    try:
        logger.info("\n" + "█" * 80)
        logger.info("STEP 1/3: Querying CELLxGENE Census for Breast Cancer Datasets")
        logger.info("█" * 80 + "\n")

        healthy_adata, tumor_adata = query_all_datasets()

        stats['healthy_cells_raw'] = healthy_adata.n_obs
        stats['tumor_cells_raw'] = tumor_adata.n_obs
        stats['genes_raw'] = healthy_adata.n_vars

        logger.info("\n✓ Step 1 completed successfully!")

    except Exception as e:
        logger.error(f"✗ Step 1 failed: {e}")
        raise

    # ========================================================================
    # Step 2: Select Highly Variable Genes (HVGs)
    # ========================================================================
    try:
        logger.info("\n" + "█" * 80)
        logger.info("STEP 2/3: Selecting Highly Variable Genes from Combined Data")
        logger.info("█" * 80 + "\n")

        healthy_filtered, tumor_filtered, hvg_list = select_hvgs(healthy_adata, tumor_adata)

        stats['healthy_cells_filtered'] = healthy_filtered.n_obs
        stats['tumor_cells_filtered'] = tumor_filtered.n_obs
        stats['n_hvgs'] = len(hvg_list)

        logger.info("\n✓ Step 2 completed successfully!")

    except Exception as e:
        logger.error(f"✗ Step 2 failed: {e}")
        raise

    # ========================================================================
    # Step 3: Map Genes to Protein Sequences
    # ========================================================================
    try:
        logger.info("\n" + "█" * 80)
        logger.info("STEP 3/3: Mapping Genes to Protein Sequences")
        logger.info("█" * 80 + "\n")

        gene_to_sequence, mapping_stats = map_genes_to_proteins(hvg_list)

        stats['n_genes_mapped'] = mapping_stats['n_mapped']
        stats['mapping_success_rate'] = mapping_stats['success_rate']
        stats['n_genes_failed'] = mapping_stats['n_total'] - mapping_stats['n_mapped']

        logger.info("\n✓ Step 3 completed successfully!")

    except Exception as e:
        logger.error(f"✗ Step 3 failed: {e}")
        raise

    # ========================================================================
    # Generate Final Report
    # ========================================================================
    elapsed_time = (datetime.now() - start_time).total_seconds()

    logger.info("\n" + "╔" + "═" * 78 + "╗")
    logger.info("║" + " " * 27 + "PIPELINE COMPLETED!" + " " * 31 + "║")
    logger.info("╚" + "═" * 78 + "╝" + "\n")

    # Print summary
    logger.info("Summary Statistics:")
    logger.info("─" * 80)
    logger.info(f"  Raw Data:")
    logger.info(f"    Healthy cells:              {stats['healthy_cells_raw']:>10,}")
    logger.info(f"    Tumor cells:                {stats['tumor_cells_raw']:>10,}")
    logger.info(f"    Protein-coding genes:       {stats['genes_raw']:>10,}")
    logger.info(f"\n  Filtered Data (HVG selected):")
    logger.info(f"    Healthy cells:              {stats['healthy_cells_filtered']:>10,}")
    logger.info(f"    Tumor cells:                {stats['tumor_cells_filtered']:>10,}")
    logger.info(f"    Highly variable genes:      {stats['n_hvgs']:>10,}")
    logger.info(f"\n  Gene-to-Protein Mapping:")
    logger.info(f"    Successfully mapped genes:  {stats['n_genes_mapped']:>10,}")
    logger.info(f"    Failed to map:              {stats['n_genes_failed']:>10,}")
    logger.info(f"    Success rate:               {stats['mapping_success_rate']:>10.1f}%")
    logger.info(f"\n  Total pipeline time:         {elapsed_time/60:>10.1f} minutes")
    logger.info("─" * 80)

    # Generate metadata summary report
    report = f"""funcCell Data Pipeline Summary
{'='*80}

Pipeline completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Total execution time: {elapsed_time/60:.1f} minutes

DATASET STATISTICS
{'─'*80}

Raw Data (after QC filtering):
  Healthy breast tissue:
    - Cells:                    {stats['healthy_cells_raw']:,}
    - Collection ID:            {config.CENSUS_CONFIG['healthy']['collection_id']}

  Breast cancer tumor:
    - Cells:                    {stats['tumor_cells_raw']:,}
    - Collection ID:            {config.CENSUS_CONFIG['tumor']['collection_id']}

  Protein-coding genes:         {stats['genes_raw']:,}

Filtered Data (HVG selection):
  Healthy cells:                {stats['healthy_cells_filtered']:,}
  Tumor cells:                  {stats['tumor_cells_filtered']:,}
  Highly variable genes:        {stats['n_hvgs']:,}

GENE-TO-PROTEIN MAPPING
{'─'*80}

  Total HVGs:                   {stats['n_hvgs']:,}
  Successfully mapped:          {stats['n_genes_mapped']:,} ({stats['mapping_success_rate']:.1f}%)
  Failed to map:                {stats['n_genes_failed']:,}

READY FOR DAY 2: ProteinBERT Embedding Generation
{'─'*80}

✓ Final gene count:             {stats['n_genes_mapped']:,} genes with protein sequences
✓ Dataset sizes known for GPU batch planning
✓ All data files saved and ready

OUTPUT FILES
{'─'*80}

Raw data:
  - {config.OUTPUT_FILES['healthy_raw']}
  - {config.OUTPUT_FILES['tumor_raw']}

Filtered data (HVG selected):
  - {config.OUTPUT_FILES['healthy_filtered']}
  - {config.OUTPUT_FILES['tumor_filtered']}

Gene-to-protein mapping:
  - {config.OUTPUT_FILES['gene_to_sequence']}
  - {config.OUTPUT_FILES['gene_list']}

Reports:
  - {config.OUTPUT_FILES['metadata_summary']}
  - {config.OUTPUT_FILES['gene_mapping_report']}

NEXT STEPS
{'─'*80}

1. Review the gene mapping report for any critical missing genes
2. Transfer data to GPU machine for ProteinBERT embedding generation
3. Use {stats['n_genes_mapped']:,} genes × batch size to plan GPU memory usage

{'='*80}
"""

    save_text(report, config.OUTPUT_FILES["metadata_summary"])
    logger.info(f"\n✓ Saved metadata summary: {config.OUTPUT_FILES['metadata_summary']}")

    logger.info("\n" + "🎉 " * 20)
    logger.info("All done! Data is ready for Day 2 ProteinBERT embedding generation.")
    logger.info("🎉 " * 20 + "\n")

    return stats


if __name__ == "__main__":
    try:
        stats = run_pipeline()
        sys.exit(0)
    except KeyboardInterrupt:
        logger.warning("\n\nPipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"\n\nPipeline failed with error: {e}")
        sys.exit(1)
