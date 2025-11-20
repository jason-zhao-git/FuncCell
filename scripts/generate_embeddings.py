#!/usr/bin/env python3
"""
Generate ProteinBERT embeddings for all mapped protein sequences.

This script:
1. Loads gene-to-sequence mapping (19,294 proteins)
2. Initializes ProteinBERT model
3. Generates 1024-dim embeddings for all sequences
4. Handles readthrough transcripts via mean pooling
5. Saves embeddings cache for downstream use

Usage:
    python scripts/generate_embeddings.py
"""

import sys
import logging
from pathlib import Path
from datetime import datetime

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from src.model.proteinbert_embeddings import ProteinBERTEmbedder, validate_embeddings
from src.preprocess_data import config as prep_config
from utils.logging_utils import setup_logger, log_step
from utils.file_utils import load_pickle, save_text

# Setup logging
logger = setup_logger(
    "embedding_generation",
    log_file=Path("data/embeddings") / f"embedding_generation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log",
    level=logging.INFO
)


@log_step(logger, "ProteinBERT Embedding Generation")
def generate_embeddings(
    seq_len: int = 1024,  # Sequences >1024 aa will be randomly subsampled
    batch_size: int = 32,
    readthrough_strategy: str = 'concat'
):
    """
    Generate ProteinBERT embeddings for all mapped genes.

    Args:
        seq_len: Maximum sequence length for encoding
        batch_size: Number of sequences to process per batch
        readthrough_strategy: Strategy for combining readthrough components
    """
    logger.info("\n")
    logger.info("╔" + "═" * 78 + "╗")
    logger.info("║" + " " * 20 + "funcCell Embedding Generation" + " " * 28 + "║")
    logger.info("║" + " " * 18 + "ProteinBERT 1024-dim Embeddings" + " " * 28 + "║")
    logger.info("╚" + "═" * 78 + "╝")
    logger.info("\n")

    start_time = datetime.now()
    stats = {}

    # ========================================================================
    # Step 1: Load Gene-to-Sequence Mapping
    # ========================================================================
    logger.info("█" * 80)
    logger.info("STEP 1/3: Loading Gene-to-Sequence Mapping")
    logger.info("█" * 80 + "\n")

    gene_to_sequence = load_pickle(prep_config.OUTPUT_FILES['gene_to_sequence'])

    stats['n_total_genes'] = len(gene_to_sequence)
    stats['n_regular_genes'] = sum(1 for v in gene_to_sequence.values() if isinstance(v, str))
    stats['n_readthrough_genes'] = sum(1 for v in gene_to_sequence.values() if isinstance(v, tuple))

    logger.info(f"Loaded {stats['n_total_genes']:,} gene sequences:")
    logger.info(f"  Regular genes:      {stats['n_regular_genes']:,}")
    logger.info(f"  Readthrough genes:  {stats['n_readthrough_genes']:,}")

    # ========================================================================
    # Step 2: Initialize ProteinBERT Model
    # ========================================================================
    logger.info("\n" + "█" * 80)
    logger.info("STEP 2/3: Initializing ProteinBERT Model")
    logger.info("█" * 80 + "\n")

    embedder = ProteinBERTEmbedder(seq_len=seq_len)

    logger.info("✓ Model initialized")
    logger.info(f"  Sequence length:  {seq_len}")
    logger.info(f"  Batch size:       {batch_size}")
    logger.info(f"  Readthrough:      {readthrough_strategy} (concatenate sequences)")

    # ========================================================================
    # Step 3: Generate Embeddings
    # ========================================================================
    logger.info("\n" + "█" * 80)
    logger.info("STEP 3/3: Generating Embeddings")
    logger.info("█" * 80 + "\n")

    cache_path = Path("data/embeddings/gene_to_embedding.pkl")
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    gene_to_embedding = embedder.generate_all_embeddings(
        gene_to_sequence,
        batch_size=batch_size,
        readthrough_strategy=readthrough_strategy,
        cache_path=cache_path,
    )

    stats['n_embeddings'] = len(gene_to_embedding)
    stats['embedding_dim'] = list(gene_to_embedding.values())[0].shape[0]

    # ========================================================================
    # Validate Embeddings
    # ========================================================================
    logger.info("\n")
    validate_embeddings(gene_to_embedding, gene_to_sequence)

    # ========================================================================
    # Generate Report
    # ========================================================================
    elapsed_time = (datetime.now() - start_time).total_seconds()
    stats['elapsed_time'] = elapsed_time

    logger.info("\n" + "╔" + "═" * 78 + "╗")
    logger.info("║" + " " * 27 + "EMBEDDING COMPLETED!" + " " * 31 + "║")
    logger.info("╚" + "═" * 78 + "╝" + "\n")

    # Print summary
    logger.info("Summary Statistics:")
    logger.info("─" * 80)
    logger.info(f"  Input Genes:")
    logger.info(f"    Total genes:                {stats['n_total_genes']:>10,}")
    logger.info(f"    Regular genes:              {stats['n_regular_genes']:>10,}")
    logger.info(f"    Readthrough genes:          {stats['n_readthrough_genes']:>10,}")
    logger.info(f"\n  Generated Embeddings:")
    logger.info(f"    Total embeddings:           {stats['n_embeddings']:>10,}")
    logger.info(f"    Embedding dimension:        {stats['embedding_dim']:>10,}")
    logger.info(f"\n  Configuration:")
    logger.info(f"    Sequence length:            {seq_len:>10,}")
    logger.info(f"    Batch size:                 {batch_size:>10,}")
    logger.info(f"    Readthrough strategy:       {readthrough_strategy:>10}")
    logger.info(f"\n  Performance:")
    logger.info(f"    Total time:                 {elapsed_time/60:>10.1f} minutes")
    logger.info(f"    Embeddings per minute:      {stats['n_embeddings']/(elapsed_time/60):>10.1f}")
    logger.info("─" * 80)

    # Generate report file
    cache_size_mb = cache_path.stat().st_size / (1024 * 1024)

    report = f"""ProteinBERT Embedding Generation Report
{'='*80}

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Execution time: {elapsed_time/60:.1f} minutes

INPUT DATA
{'─'*80}

  Gene-to-sequence mapping:    {prep_config.OUTPUT_FILES['gene_to_sequence']}
  Total genes:                  {stats['n_total_genes']:,}
    - Regular genes:            {stats['n_regular_genes']:,}
    - Readthrough genes:        {stats['n_readthrough_genes']:,}

PROTEINBERT CONFIGURATION
{'─'*80}

  Model:                        ProteinBERT (pretrained on 106M proteins)
  Sequence length:              {seq_len} tokens
  Batch size:                   {batch_size}
  Readthrough strategy:         {readthrough_strategy} pooling
  Embedding dimension:          {stats['embedding_dim']}

RESULTS
{'─'*80}

  Total embeddings generated:   {stats['n_embeddings']:,}
  Coverage:                     {100 * stats['n_embeddings'] / stats['n_total_genes']:.1f}%
  All embeddings validated:     ✓

  Cache file:                   {cache_path}
  Cache size:                   {cache_size_mb:.1f} MB

PERFORMANCE
{'─'*80}

  Total runtime:                {elapsed_time/60:.1f} minutes
  Embeddings per minute:        {stats['n_embeddings']/(elapsed_time/60):.1f}
  Average time per embedding:   {elapsed_time/stats['n_embeddings']:.3f} seconds

NEXT STEPS
{'─'*80}

1. Load embeddings cache: gene_to_embedding = pickle.load(open('{cache_path}', 'rb'))
2. Combine with single-cell expression data for cell-level representations
3. Train cancer vs. healthy classifier

{'='*80}
"""

    report_path = Path("data/embeddings/embedding_generation_report.txt")
    save_text(report, report_path)
    logger.info(f"\n✓ Saved report: {report_path}")

    logger.info("\n" + "🎉 " * 20)
    logger.info("Embeddings successfully generated!")
    logger.info("Ready for downstream cancer classification.")
    logger.info("🎉 " * 20 + "\n")

    return gene_to_embedding


if __name__ == "__main__":
    try:
        gene_to_embedding = generate_embeddings(
            seq_len=1024,  # Sequences >1024 aa will be randomly subsampled
            batch_size=32,
            readthrough_strategy='concat'
        )
        sys.exit(0)
    except KeyboardInterrupt:
        logger.warning("\n\nEmbedding generation interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"\n\nEmbedding generation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
