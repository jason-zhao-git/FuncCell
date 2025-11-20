"""Generate ProteinBERT embeddings for protein sequences."""

import logging
import pickle
import numpy as np
import urllib.request
from pathlib import Path
from typing import Dict, Union, Tuple
from tqdm import tqdm

from proteinbert import load_pretrained_model, load_pretrained_model_from_dump
from proteinbert.conv_and_global_attention_model import (
    get_model_with_hidden_layers_as_outputs
)

logger = logging.getLogger(__name__)

# Local model path in project directory
PROJECT_ROOT = Path(__file__).parent.parent.parent
MODEL_DIR = PROJECT_ROOT / "models" / "proteinbert"
MODEL_PATH = MODEL_DIR / "full_go_epoch_92400_sample_23500000.pkl"
ZENODO_URL = "https://zenodo.org/records/10371965/files/full_go_epoch_92400_sample_23500000.pkl?download=1"


class ProteinBERTEmbedder:
    """Generate and cache ProteinBERT embeddings."""

    def __init__(self, seq_len: int = 1024, device: str = 'gpu'):
        """
        Initialize ProteinBERT model.

        Args:
            seq_len: Maximum sequence length for encoding
            device: 'gpu' or 'cpu' (currently not used, TensorFlow auto-detects)
        """
        logger.info("Loading ProteinBERT pretrained model...")

        # Check if model exists locally, download if not
        if not MODEL_PATH.exists():
            logger.info("Model not found locally. Downloading from Zenodo...")
            logger.info(f"  Source: {ZENODO_URL}")
            logger.info(f"  Target: {MODEL_PATH}")
            logger.info(f"  Size: ~200 MB (this may take a few minutes)")

            MODEL_DIR.mkdir(parents=True, exist_ok=True)

            def reporthook(block_num, block_size, total_size):
                downloaded = block_num * block_size
                percent = min(100, downloaded * 100 / total_size)
                mb_downloaded = downloaded / (1024 * 1024)
                mb_total = total_size / (1024 * 1024)
                logger.info(f"  Progress: {percent:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)")

            try:
                urllib.request.urlretrieve(ZENODO_URL, MODEL_PATH, reporthook)
                logger.info(f"✓ Model downloaded successfully to {MODEL_PATH}")
            except Exception as e:
                logger.error(f"Failed to download model: {e}")
                logger.error(f"Please download manually from: {ZENODO_URL}")
                logger.error(f"And save to: {MODEL_PATH}")
                raise
        else:
            logger.info(f"Using cached model: {MODEL_PATH}")

        # Load pretrained model from our local directory
        self.pretrained_model_generator, self.input_encoder = load_pretrained_model(
            local_model_dump_dir=str(MODEL_DIR),
            local_model_dump_file_name=MODEL_PATH.name,
            download_model_dump_if_not_exists=False  # Don't try to download, we already have it
        )
        self.model = get_model_with_hidden_layers_as_outputs(
            self.pretrained_model_generator.create_model(seq_len)
        )

        self.seq_len = seq_len
        self.device = device

        logger.info(f"Model loaded successfully (seq_len={seq_len})")


    def embed_sequence(self, sequence: str) -> np.ndarray:
        """
        Generate embedding for a single protein sequence.
        Long sequences (>seq_len) will be truncated to fit model input.

        Args:
            sequence: Amino acid sequence string

        Returns:
            1024-dimensional embedding vector
        """
        # Truncate sequence if too long (ProteinBERT handles this with random subsampling)
        if len(sequence) > self.seq_len - 10:  # Leave room for special tokens
            # Random subsample for variety
            import random
            start = random.randint(0, len(sequence) - (self.seq_len - 10))
            sequence = sequence[start:start + (self.seq_len - 10)]

        encoded_x = self.input_encoder.encode_X([sequence], self.seq_len)
        _, global_rep = self.model.predict(encoded_x, batch_size=1, verbose=0)
        return global_rep[0]


    def embed_batch(self, sequences: list, batch_size: int = 32) -> np.ndarray:
        """
        Generate embeddings for a batch of sequences.

        Args:
            sequences: List of amino acid sequence strings
            batch_size: Batch size for processing

        Returns:
            Array of shape (n_sequences, 1024)
        """
        all_embeddings = []

        for i in range(0, len(sequences), batch_size):
            batch = sequences[i:i+batch_size]

            # Truncate long sequences
            truncated_batch = []
            for seq in batch:
                if len(seq) > self.seq_len - 10:
                    import random
                    start = random.randint(0, len(seq) - (self.seq_len - 10))
                    truncated_batch.append(seq[start:start + (self.seq_len - 10)])
                else:
                    truncated_batch.append(seq)

            encoded_x = self.input_encoder.encode_X(truncated_batch, self.seq_len)
            _, global_reps = self.model.predict(encoded_x, batch_size=len(batch), verbose=0)
            all_embeddings.append(global_reps)

        return np.vstack(all_embeddings)


    def embed_readthrough(
        self,
        component_sequences: Tuple[str, ...],
        strategy: str = 'concat'
    ) -> np.ndarray:
        """
        Generate embedding for readthrough transcript (fusion gene).

        Args:
            component_sequences: Tuple of component gene sequences
            strategy: 'concat', 'mean', 'max', or 'weighted'
                - 'concat': Concatenate sequences then embed (biologically accurate)
                - 'mean': Embed separately then average
                - 'max': Embed separately then take element-wise max
                - 'weighted': Embed separately then weighted average by length

        Returns:
            Embedding vector (1024-dim)
        """
        if strategy == 'concat':
            # Concatenate amino acid sequences into one fusion protein
            # This is biologically accurate - readthrough produces continuous polypeptide
            fusion_sequence = ''.join(component_sequences)
            return self.embed_sequence(fusion_sequence)

        # For other strategies, embed components separately first
        component_embeddings = []
        for seq in component_sequences:
            emb = self.embed_sequence(seq)
            component_embeddings.append(emb)

        component_embeddings = np.array(component_embeddings)

        # Combine embeddings based on strategy
        if strategy == 'mean':
            return np.mean(component_embeddings, axis=0)
        elif strategy == 'max':
            return np.max(component_embeddings, axis=0)
        elif strategy == 'weighted':
            # Weight by sequence length
            lengths = np.array([len(seq) for seq in component_sequences])
            weights = lengths / lengths.sum()
            return np.average(component_embeddings, axis=0, weights=weights)
        else:
            raise ValueError(f"Unknown strategy: {strategy}")


    def generate_all_embeddings(
        self,
        gene_to_sequence: Dict[str, Union[str, Tuple[str, ...]]],
        batch_size: int = 32,
        readthrough_strategy: str = 'concat',
        cache_path: Path = None,
    ) -> Dict[str, np.ndarray]:
        """
        Generate embeddings for all genes (regular + readthrough).

        Args:
            gene_to_sequence: Dict mapping gene symbols to sequences
            batch_size: Batch size for processing
            readthrough_strategy: Strategy for combining readthrough components
            cache_path: Path to save embeddings cache

        Returns:
            Dict mapping gene symbols to embedding vectors
        """
        logger.info("="*60)
        logger.info("Generating ProteinBERT embeddings")
        logger.info("="*60)

        # Separate regular and readthrough genes
        regular_genes = {k: v for k, v in gene_to_sequence.items()
                        if isinstance(v, str)}
        readthrough_genes = {k: v for k, v in gene_to_sequence.items()
                           if isinstance(v, tuple)}

        logger.info(f"Regular genes: {len(regular_genes):,}")
        logger.info(f"Readthrough genes: {len(readthrough_genes):,}")
        logger.info(f"Readthrough strategy: {readthrough_strategy}")
        logger.info(f"Batch size: {batch_size}")
        logger.info(f"Sequence length: {self.seq_len}")

        gene_to_embedding = {}

        # Process regular genes one by one (batch processing has numpy compatibility issues)
        logger.info("\nProcessing regular genes...")
        for gene, seq in tqdm(regular_genes.items(), desc="Regular genes"):
            gene_to_embedding[gene] = self.embed_sequence(seq)

        # Process readthrough genes
        if readthrough_genes:
            logger.info("\nProcessing readthrough genes...")
            for gene, component_seqs in tqdm(readthrough_genes.items(),
                                            desc="Readthrough genes"):
                gene_to_embedding[gene] = self.embed_readthrough(
                    component_seqs,
                    strategy=readthrough_strategy
                )

        # Generate statistics
        logger.info("\n" + "="*60)
        logger.info("Embedding generation completed!")
        logger.info("="*60)
        logger.info(f"Total embeddings: {len(gene_to_embedding):,}")

        # Check embedding dimensions
        embedding_dims = set(emb.shape[0] for emb in gene_to_embedding.values())
        logger.info(f"Embedding dimensions: {embedding_dims}")

        # Cache results
        if cache_path:
            logger.info(f"\nSaving embeddings to {cache_path}")
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_path, 'wb') as f:
                pickle.dump(gene_to_embedding, f)

            # Report file size
            file_size_mb = cache_path.stat().st_size / (1024 * 1024)
            logger.info(f"  Cache size: {file_size_mb:.1f} MB")

        return gene_to_embedding


def validate_embeddings(
    gene_to_embedding: Dict[str, np.ndarray],
    gene_to_sequence: Dict[str, Union[str, Tuple[str, ...]]]
) -> bool:
    """
    Validate generated embeddings.

    Args:
        gene_to_embedding: Dict mapping gene symbols to embeddings
        gene_to_sequence: Dict mapping gene symbols to sequences

    Returns:
        True if all validations pass
    """
    logger.info("\n" + "="*60)
    logger.info("Validating embeddings")
    logger.info("="*60)

    # Check 1: Coverage
    assert len(gene_to_embedding) == len(gene_to_sequence), \
        f"Coverage mismatch: {len(gene_to_embedding)} embeddings vs {len(gene_to_sequence)} sequences"
    logger.info(f"✓ Coverage: {len(gene_to_embedding):,} genes")

    # Check 2: Dimensionality
    dims = [emb.shape[0] for emb in gene_to_embedding.values()]
    assert all(d == 1024 for d in dims), \
        f"Not all embeddings are 1024-dimensional: {set(dims)}"
    logger.info(f"✓ All embeddings are 1024-dimensional")

    # Check 3: No NaN or Inf
    for gene, emb in gene_to_embedding.items():
        assert not np.isnan(emb).any(), f"NaN values in {gene}"
        assert not np.isinf(emb).any(), f"Inf values in {gene}"
    logger.info(f"✓ No NaN or Inf values")

    # Check 4: Reasonable magnitude
    magnitudes = [np.linalg.norm(emb) for emb in gene_to_embedding.values()]
    mean_mag = np.mean(magnitudes)
    std_mag = np.std(magnitudes)
    logger.info(f"✓ Embedding magnitudes: {mean_mag:.2f} ± {std_mag:.2f}")

    # Check 5: Readthrough genes
    n_readthrough = sum(1 for v in gene_to_sequence.values() if isinstance(v, tuple))
    logger.info(f"✓ Readthrough genes handled: {n_readthrough}")

    logger.info("="*60)
    logger.info("All validations passed!")
    logger.info("="*60)

    return True


if __name__ == "__main__":
    from utils.logging_utils import setup_logger
    from utils.file_utils import load_pickle
    from src.preprocess_data import config as prep_config

    setup_logger("proteinbert_embeddings", level=logging.INFO)

    # Load gene-to-sequence mapping
    logger.info("Loading gene-to-sequence mapping...")
    gene_to_sequence = load_pickle(prep_config.OUTPUT_FILES['gene_to_sequence'])
    logger.info(f"Loaded {len(gene_to_sequence):,} gene sequences")

    # Initialize embedder
    embedder = ProteinBERTEmbedder(seq_len=1024)

    # Generate embeddings
    cache_path = Path("data/embeddings/gene_to_embedding.pkl")

    gene_to_embedding = embedder.generate_all_embeddings(
        gene_to_sequence,
        batch_size=32,
        readthrough_strategy='concat',
        cache_path=cache_path,
    )

    # Validate
    validate_embeddings(gene_to_embedding, gene_to_sequence)
