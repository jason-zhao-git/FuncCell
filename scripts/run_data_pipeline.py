#!/usr/bin/env python3
"""
Main script to run the Day 1 data preprocessing pipeline.

This script:
1. Loads all protein-coding genes from BioMart
2. Maps genes to protein sequences via MyGene.info and UniProt
3. Generates comprehensive reports

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
from src.preprocess_data.gene_list import get_protein_coding_genes
from src.preprocess_data.gene_mapping import map_genes_to_proteins
from utils.logging_utils import setup_logger, log_step
from utils.file_utils import save_text, save_list

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
    logger.info("║" + " " * 15 + "Protein-Coding Gene to Sequence Mapping" + " " * 24 + "║")
    logger.info("╚" + "═" * 78 + "╝")
    logger.info("\n")

    # Track overall statistics
    stats = {}
    start_time = datetime.now()

    # ========================================================================
    # Step 1: Load Protein-Coding Genes from BioMart
    # ========================================================================
    try:
        logger.info("\n" + "█" * 80)
        logger.info("STEP 1/2: Loading Protein-Coding Genes from BioMart")
        logger.info("█" * 80 + "\n")

        gene_list = get_protein_coding_genes()

        stats['n_genes_total'] = len(gene_list)

        # Save gene list
        save_list(gene_list, config.OUTPUT_FILES["gene_list"])
        logger.info(f"Saved {len(gene_list):,} genes to: {config.OUTPUT_FILES['gene_list']}")

        logger.info("\n✓ Step 1 completed successfully!")

    except Exception as e:
        logger.error(f"✗ Step 1 failed: {e}")
        raise

    # ========================================================================
    # Step 2: Map Genes to Protein Sequences
    # ========================================================================
    try:
        logger.info("\n" + "█" * 80)
        logger.info("STEP 2/2: Mapping Genes to Protein Sequences")
        logger.info("█" * 80 + "\n")

        gene_to_sequence, mapping_stats = map_genes_to_proteins(gene_list)

        stats['n_genes_mapped'] = mapping_stats['n_mapped']
        stats['mapping_success_rate'] = mapping_stats['success_rate']
        stats['n_genes_failed'] = mapping_stats['n_total'] - mapping_stats['n_mapped']

        logger.info("\n✓ Step 2 completed successfully!")

    except Exception as e:
        logger.error(f"✗ Step 2 failed: {e}")
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
    logger.info(f"  Protein-Coding Genes (BioMart):")
    logger.info(f"    Total genes:                {stats['n_genes_total']:>10,}")
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

GENE SOURCE
{'─'*80}

  BioMart protein-coding genes: {stats['n_genes_total']:,}
  Source file:                  {config.GENE_FILTER_PARAMS['biomart_file']}

GENE-TO-PROTEIN MAPPING
{'─'*80}

  Total genes:                  {stats['n_genes_total']:,}
  Successfully mapped:          {stats['n_genes_mapped']:,} ({stats['mapping_success_rate']:.1f}%)
  Failed to map:                {stats['n_genes_failed']:,}

READY FOR DAY 2: ProteinBERT Embedding Generation
{'─'*80}

✓ Final gene count:             {stats['n_genes_mapped']:,} genes with protein sequences
✓ Universal embedding cache (works with any dataset)
✓ Ready to generate ProteinBERT embeddings on GPU

OUTPUT FILES
{'─'*80}

Gene-to-protein mapping:
  - {config.OUTPUT_FILES['gene_to_sequence']}
  - {config.OUTPUT_FILES['gene_list']}

Reports:
  - {config.OUTPUT_FILES['metadata_summary']}
  - {config.OUTPUT_FILES['gene_mapping_report']}

NEXT STEPS
{'─'*80}

1. Review the gene mapping report for any critical missing genes
2. Transfer gene_to_sequence.pkl to GPU machine
3. Generate ProteinBERT embeddings for {stats['n_genes_mapped']:,} protein sequences
4. When querying CELLxGENE Census, filter to genes in the embedding cache

{'='*80}
"""

    save_text(report, config.OUTPUT_FILES["metadata_summary"])
    logger.info(f"\n✓ Saved metadata summary: {config.OUTPUT_FILES['metadata_summary']}")

    logger.info("\n" + "🎉 " * 20)
    logger.info("All done! Gene-to-protein mapping complete.")
    logger.info("Ready for ProteinBERT embedding generation.")
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
