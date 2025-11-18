"""Map genes to protein sequences using MyGene and UniProt APIs."""

import logging
import time
from typing import Dict, List, Tuple
import requests
import mygene
from tqdm import tqdm

from . import config
from utils.file_utils import save_pickle, save_list, save_text

logger = logging.getLogger(__name__)


def query_mygene(gene_symbols: List[str], batch_size: int = 100) -> Dict[str, str]:
    """
    Query MyGene.info to map gene symbols to UniProt IDs.

    Args:
        gene_symbols: List of gene symbols
        batch_size: Number of genes per batch request

    Returns:
        Dictionary mapping gene symbol to UniProt ID
    """
    logger.info(f"Querying MyGene.info for {len(gene_symbols)} genes...")

    mg = mygene.MyGeneInfo()
    gene_to_uniprot = {}

    # Process in batches
    for i in tqdm(range(0, len(gene_symbols), batch_size), desc="Querying MyGene"):
        batch = gene_symbols[i:i + batch_size]

        try:
            # Query MyGene
            results = mg.querymany(
                batch,
                scopes='symbol',
                fields='uniprot.Swiss-Prot',
                species='human',
                returnall=True,
            )

            # Parse results
            for result in results['out']:
                gene_symbol = result.get('query')
                if 'uniprot' in result and 'Swiss-Prot' in result['uniprot']:
                    # Get first Swiss-Prot ID (canonical)
                    uniprot_id = result['uniprot']['Swiss-Prot']
                    if isinstance(uniprot_id, list):
                        uniprot_id = uniprot_id[0]
                    gene_to_uniprot[gene_symbol] = uniprot_id

            # Be nice to the API
            time.sleep(0.2)

        except Exception as e:
            logger.warning(f"Error querying batch {i//batch_size + 1}: {e}")
            continue

    logger.info(f"Mapped {len(gene_to_uniprot)}/{len(gene_symbols)} genes to UniProt IDs")
    logger.info(f"Mapping success rate: {100 * len(gene_to_uniprot) / len(gene_symbols):.1f}%")

    return gene_to_uniprot


def fetch_uniprot_sequences(
    uniprot_ids: List[str],
    batch_size: int = 100,
    max_retries: int = 3,
) -> Dict[str, str]:
    """
    Fetch protein sequences from UniProt REST API.

    Args:
        uniprot_ids: List of UniProt IDs
        batch_size: Number of IDs per batch request
        max_retries: Maximum retry attempts for failed requests

    Returns:
        Dictionary mapping UniProt ID to protein sequence
    """
    logger.info(f"Fetching protein sequences from UniProt for {len(uniprot_ids)} IDs...")

    uniprot_to_seq = {}
    base_url = "https://rest.uniprot.org/uniprotkb/stream"

    # Process in batches
    for i in tqdm(range(0, len(uniprot_ids), batch_size), desc="Fetching UniProt"):
        batch = uniprot_ids[i:i + batch_size]
        query = " OR ".join([f"accession:{uid}" for uid in batch])

        # Retry logic
        for attempt in range(max_retries):
            try:
                params = {
                    'query': query,
                    'format': 'fasta',
                    'compressed': 'false',
                }

                response = requests.get(
                    base_url,
                    params=params,
                    timeout=config.GENE_MAPPING_PARAMS["timeout"],
                )
                response.raise_for_status()

                # Parse FASTA response
                fasta_text = response.text
                current_id = None
                current_seq = []

                for line in fasta_text.split('\n'):
                    if line.startswith('>'):
                        # Save previous sequence
                        if current_id and current_seq:
                            uniprot_to_seq[current_id] = ''.join(current_seq)

                        # Parse header: >sp|P12345|GENE_HUMAN ...
                        parts = line.split('|')
                        if len(parts) >= 2:
                            current_id = parts[1]
                        current_seq = []
                    elif line.strip():
                        current_seq.append(line.strip())

                # Save last sequence
                if current_id and current_seq:
                    uniprot_to_seq[current_id] = ''.join(current_seq)

                # Success, break retry loop
                break

            except Exception as e:
                if attempt < max_retries - 1:
                    logger.warning(f"Retry {attempt + 1}/{max_retries} for batch {i//batch_size + 1}: {e}")
                    time.sleep(2 ** attempt)  # Exponential backoff
                else:
                    logger.error(f"Failed to fetch batch {i//batch_size + 1} after {max_retries} attempts")

        # Be nice to the API
        time.sleep(0.3)

    logger.info(f"Fetched {len(uniprot_to_seq)}/{len(uniprot_ids)} protein sequences")
    logger.info(f"Fetch success rate: {100 * len(uniprot_to_seq) / len(uniprot_ids):.1f}%")

    return uniprot_to_seq


def map_genes_to_proteins(gene_list: List[str]) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Map gene symbols to protein sequences.

    Args:
        gene_list: List of gene symbols

    Returns:
        Tuple of (gene_to_sequence dict, mapping_stats dict)
    """
    logger.info("="*60)
    logger.info("Mapping genes to protein sequences")
    logger.info("="*60)

    # Step 1: Gene symbol → UniProt ID
    gene_to_uniprot = query_mygene(
        gene_list,
        batch_size=config.GENE_MAPPING_PARAMS["batch_size"]
    )

    # Step 2: UniProt ID → Protein sequence
    uniprot_ids = list(gene_to_uniprot.values())
    uniprot_to_seq = fetch_uniprot_sequences(
        uniprot_ids,
        batch_size=config.GENE_MAPPING_PARAMS["batch_size"],
        max_retries=config.GENE_MAPPING_PARAMS["max_retries"],
    )

    # Step 3: Combine mappings: Gene → Sequence
    gene_to_sequence = {}
    for gene, uniprot_id in gene_to_uniprot.items():
        if uniprot_id in uniprot_to_seq:
            gene_to_sequence[gene] = uniprot_to_seq[uniprot_id]

    # Calculate statistics
    n_total = len(gene_list)
    n_mapped = len(gene_to_sequence)
    success_rate = 100 * n_mapped / n_total

    logger.info("\n" + "="*60)
    logger.info("Gene mapping completed!")
    logger.info("="*60)
    logger.info(f"Total genes:              {n_total:,}")
    logger.info(f"Successfully mapped:      {n_mapped:,}")
    logger.info(f"Failed to map:            {n_total - n_mapped:,}")
    logger.info(f"Overall success rate:     {success_rate:.1f}%")

    # Find missing genes
    missing_genes = [g for g in gene_list if g not in gene_to_sequence]

    # Calculate sequence length statistics
    seq_lengths = [len(seq) for seq in gene_to_sequence.values()]
    if seq_lengths:
        logger.info(f"\nProtein sequence statistics:")
        logger.info(f"  Mean length:   {sum(seq_lengths) / len(seq_lengths):.0f} aa")
        logger.info(f"  Median length: {sorted(seq_lengths)[len(seq_lengths)//2]:.0f} aa")
        logger.info(f"  Min length:    {min(seq_lengths)} aa")
        logger.info(f"  Max length:    {max(seq_lengths)} aa")

    # Save results
    logger.info("\nSaving results...")
    save_pickle(gene_to_sequence, config.OUTPUT_FILES["gene_to_sequence"])

    # Save final gene list (only successfully mapped genes)
    mapped_genes = sorted(gene_to_sequence.keys())
    save_list(mapped_genes, config.OUTPUT_FILES["gene_list"])

    # Generate mapping report
    report = f"""Gene-to-Protein Mapping Report
{'='*60}

Summary:
  Total input genes:        {n_total:,}
  Successfully mapped:      {n_mapped:,} ({success_rate:.1f}%)
  Failed to map:            {n_total - n_mapped:,}

Protein Sequence Statistics:
  Mean length:              {sum(seq_lengths) / len(seq_lengths):.0f} aa
  Median length:            {sorted(seq_lengths)[len(seq_lengths)//2]:.0f} aa
  Range:                    {min(seq_lengths)} - {max(seq_lengths)} aa

Missing Genes ({len(missing_genes)}):
{chr(10).join('  - ' + g for g in missing_genes[:50])}
{'  ... and ' + str(len(missing_genes) - 50) + ' more' if len(missing_genes) > 50 else ''}

Output Files:
  - Gene-to-sequence mapping:  {config.OUTPUT_FILES['gene_to_sequence']}
  - Final gene list:           {config.OUTPUT_FILES['gene_list']}
  - This report:               {config.OUTPUT_FILES['gene_mapping_report']}
"""

    save_text(report, config.OUTPUT_FILES["gene_mapping_report"])
    logger.info(f"  Saved mapping report: {config.OUTPUT_FILES['gene_mapping_report']}")

    mapping_stats = {
        'n_total': n_total,
        'n_mapped': n_mapped,
        'success_rate': success_rate,
        'missing_genes': missing_genes,
    }

    return gene_to_sequence, mapping_stats


if __name__ == "__main__":
    from utils.logging_utils import setup_logger
    from utils.file_utils import load_list

    setup_logger("gene_mapping", level=logging.INFO)

    # Load gene list from HVG selection
    gene_list = load_list(config.OUTPUT_FILES["gene_list"])
    logger.info(f"Loaded {len(gene_list)} genes from HVG selection")

    # Map genes to proteins
    map_genes_to_proteins(gene_list)
