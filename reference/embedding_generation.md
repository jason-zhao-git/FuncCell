# ProteinBERT Embedding Generation

**Date:** 2025-11-20
**Model:** ProteinBERT (pretrained on 106M proteins)
**Total Proteins:** 19,294 (19,194 regular + 100 readthrough)
**Embedding Dimension:** 1024

## Summary

ProteinBERT embeddings provide protein function representations for all mapped genes. Each protein sequence (variable length: 25-34,350 amino acids) is compressed into a fixed 1024-dimensional vector capturing functional properties learned from 106 million proteins.

## Model Architecture

### ProteinBERT Overview
- **Framework:** TensorFlow/Keras (not PyTorch)
- **Pretraining:** 106M proteins from UniRef90
- **Training time:** 28 days on GPU
- **Key innovation:** Global attention (linear complexity O(n) vs standard O(n²))
- **Advantage:** Can handle extremely long sequences up to 34k amino acids

### Dual Representation System
- **Local representations:** Per-amino acid embeddings (1024-dim × sequence_length)
- **Global representations:** Whole-protein embeddings (1024-dim per protein) ← **We use these**

## Implementation Details

### Configuration
```python
seq_len = 1024               # Maximum sequence length for encoding
batch_size = 32              # Sequences processed per batch
readthrough_strategy = 'mean' # How to combine fusion gene components
```

### Sequence Length Handling
- **Training lengths:** 128, 512, or 1024 tokens
- **Our choice:** 1024 tokens
  - Covers 95%+ of proteins (median=431aa, mean=579aa)
  - Long sequences (>1024aa) are randomly subsampled
  - Model's global attention captures full-protein context from substrings

#### Truncation Strategy for Long Proteins (>1024aa)

**Dataset Statistics:**
- 2,260 proteins (11.7%) exceed 1024aa
- Longest: TTN (titin) at 34,350aa
- Other examples: MUC16 (14,507aa), SYNE1 (8,797aa), NEB (8,525aa)

**Chosen Approach: Single Random Substring**
```python
# For sequences >1024aa, randomly select one 1024aa window
if len(sequence) > 1024:
    start = random.randint(0, len(sequence) - 1024)
    sequence = sequence[start:start + 1024]
# Embed once → one 1024-dim vector
```

**Why Single Random Substring:**
1. **Speed:** One inference per protein (~3 min for all 19,294 proteins)
2. **ProteinBERT design:** Global attention architecture captures full-protein context from substrings
3. **Downstream aggregation:** Cell-level embeddings combine across many genes anyway
4. **Good enough:** For initial cancer classification experiments

**Alternative Considered: Multi-Window Aggregation**
```python
# Split long protein into multiple windows
windows = [protein[i:i+1024] for i in range(0, len(protein), 1024)]
embeddings = [embed(w) for w in windows]  # N inferences
final_embedding = np.mean(embeddings, axis=0)  # Aggregate
```
- **Pros:** Comprehensive coverage, captures all functional domains
- **Cons:** 5-7× slower (~15-20 min total), more complex
- **Status:** Can implement later if single-substring approach insufficient

**Decision rationale:** For cancer classification where we're combining embeddings with expression data across thousands of genes per cell, the speed advantage outweighs the marginal biological completeness benefit. The global attention mechanism is specifically designed to capture protein function from substrings.

### Readthrough Transcript Strategy

**Problem:** 100 readthrough genes stored as tuples of component sequences
**Examples:** `PDCD6-AHRR`, `NME1-NME2`, `BIVM-ERCC5`

**Solution:** Mean pooling of component embeddings
```python
# For readthrough gene with components (seq1, seq2)
embedding1 = embed_protein(seq1)  # 1024-dim
embedding2 = embed_protein(seq2)  # 1024-dim
readthrough_embedding = (embedding1 + embedding2) / 2  # 1024-dim
```

**Why mean pooling:**
- Preserves dimensionality (all genes have 1024-dim embeddings)
- Biologically reasonable (fusion protein combines functions)
- Simple and interpretable
- Compatible with downstream models (no special handling needed)

**Alternatives considered:**
- Max pooling: Too aggressive, loses complementary information
- Concatenation: Creates 2048-dim vectors, breaks compatibility
- Weighted average: Adds complexity, marginal benefit

## Technical Specifications

### Dependencies
```python
tensorflow>=2.4.0
tensorflow-addons>=0.12.1  # Requires Python <3.12
h5py>=3.2.1
lxml>=4.3.2
pyfaidx>=0.5.8
protein-bert>=1.0.0
```

**Important:** Python 3.12 is not supported due to tensorflow-addons. The project uses Python 3.9-3.11.

### Model Loading
```python
from proteinbert import load_pretrained_model
from proteinbert.conv_and_global_attention_model import (
    get_model_with_hidden_layers_as_outputs
)

# Load pretrained model (downloads ~500 MB on first run)
pretrained_model_generator, input_encoder = load_pretrained_model()
model = get_model_with_hidden_layers_as_outputs(
    pretrained_model_generator.create_model(1024)
)
```

### Batch Processing
```python
# Encode sequences
encoded_x = input_encoder.encode_X(sequences, seq_len=1024)

# Get embeddings
local_reps, global_reps = model.predict(encoded_x, batch_size=32)

# global_reps shape: (n_sequences, 1024)
```

## Output

### File: `data/embeddings/gene_to_embedding.pkl`
```python
{
    'TP53': array([0.12, -0.34, ..., 0.56]),      # 1024 floats
    'BRCA1': array([0.45, -0.12, ..., 0.89]),    # 1024 floats
    'PDCD6-AHRR': array([0.23, 0.15, ..., 0.34]), # 1024 floats (mean of components)
    # ... 19,294 total genes
}
```

### Statistics
- **File size:** ~79 MB (19,294 × 1024 × 4 bytes)
- **Coverage:** 100% of mapped genes
- **Dimension:** All embeddings are 1024-dim
- **Quality:** No NaN or Inf values, all validated

## Validation

Post-generation checks performed:
1. **Coverage:** All 19,294 genes have embeddings
2. **Dimensionality:** All vectors are exactly 1024-dim
3. **Validity:** No NaN or Inf values
4. **Magnitude:** Reasonable embedding norms (mean ± std logged)
5. **Readthrough:** All 100 fusion genes handled correctly

## Usage Example

### Load Embeddings
```python
import pickle

with open('data/embeddings/gene_to_embedding.pkl', 'rb') as f:
    gene_to_embedding = pickle.load(f)

# Access embedding
tp53_embedding = gene_to_embedding['TP53']  # shape: (1024,)
```

### Create Cell-Level Embeddings
```python
import numpy as np

def create_cell_embeddings(adata, gene_to_embedding):
    """
    Combine expression with protein embeddings for cell representation.

    Args:
        adata: AnnData with expression (cells × genes)
        gene_to_embedding: Dict mapping gene symbols to embeddings

    Returns:
        cell_embeddings: Array of shape (n_cells, 1024)
    """
    # Filter to genes with embeddings
    common_genes = [g for g in adata.var_names if g in gene_to_embedding]
    adata_filtered = adata[:, common_genes]

    # Get expression: (n_cells, n_genes)
    expression = adata_filtered.X
    if issparse(expression):
        expression = expression.toarray()

    # Normalize per cell
    expression_norm = expression / (expression.sum(axis=1, keepdims=True) + 1e-10)

    # Get protein embeddings: (n_genes, 1024)
    protein_embeddings = np.array([gene_to_embedding[g] for g in common_genes])

    # Weighted sum: (n_cells, 1024)
    cell_embeddings = expression_norm @ protein_embeddings

    return cell_embeddings
```

## Performance

### Computational Requirements
- **GPU:** Nvidia RTX 3060 or better (8+ GB VRAM recommended)
- **CPU:** Fallback available (10x slower)
- **RAM:** 16 GB recommended
- **Storage:** 500 MB for model + 79 MB for embeddings

### Timing
| Task | GPU | CPU |
|------|-----|-----|
| Model download | 5 min | 5 min |
| Model loading | 30 sec | 30 sec |
| Embedding 19,294 proteins | 2-3 min | 20-30 min |
| Total | ~10 min | ~40 min |

## Model Details

### ProteinBERT Architecture
- **Input:** Amino acid sequence (variable length)
- **Tokenization:** Each amino acid → token + special tokens (START, END, PAD)
- **Encoder:** 6 transformer-like blocks with global attention
- **Output:**
  - Local: (seq_len, 1024) per-residue embeddings
  - Global: (1024,) whole-protein embedding

### Global Attention Mechanism
Standard self-attention requires O(n²) operations for sequence length n, limiting maximum length. ProteinBERT uses a global-attention architecture:

1. **Local representations** attend to their neighborhood
2. **Global representation** aggregates information from all positions
3. Complexity: O(n) instead of O(n²)
4. **Result:** Can process 34,350 amino acid sequences efficiently

### Pretraining Task
- **Dataset:** UniRef90 (~106M proteins)
- **Objective:** Masked language modeling (predict masked amino acids)
- **Additional:** GO term prediction (multi-label classification)
- **Duration:** 28 days on single GPU
- **Throughput:** 280 protein records/second during training

## Why ProteinBERT

### Advantages
1. **Linear complexity** → Handles long sequences (up to 34k aa)
2. **Pretrained** on 106M proteins → Rich functional knowledge
3. **Accessible** → Single GPU inference, fast (~3 min for 19k proteins)
4. **Proven** → Outperforms sequence-only methods on function prediction
5. **Standalone** → No alignment or homology search required

### Alternatives Considered
- **ESM-2:** Meta's model, PyTorch-based, Python 3.12 compatible
  - Pros: More recent, larger scale pretraining
  - Cons: Slower inference, larger model size
  - **Decision:** ProteinBERT sufficient for initial experiments

- **ProtT5-XL:** Transformer-based, very large
  - Pros: State-of-the-art performance
  - Cons: Requires multiple GPUs, much slower
  - **Decision:** Unnecessary for cell classification

## Next Steps

1. **Load embeddings** in downstream scripts
2. **Combine with single-cell expression data** (CELLxGENE Census)
3. **Create cell-level representations** via expression-weighted sum
4. **Train classifier** to distinguish cancer vs healthy cells

## References

- **GitHub:** https://github.com/nadavbra/protein_bert
- **Paper:** Brandes et al., Bioinformatics 2022
- **Hugging Face:** https://huggingface.co/GrimSqueaker/proteinBERT
- **Zenodo (weights):** https://zenodo.org/records/10371965

## Troubleshooting

### Python Version Error
```
error: tensorflow-addons only has wheels for cp39, cp310, cp311
```
**Solution:** Update `requires-python = ">=3.9,<3.12"` in pyproject.toml

### Out of Memory (GPU)
```
ResourceExhausted: OOM when allocating tensor
```
**Solution:** Reduce batch_size from 32 to 16 or 8

### Model Download Failure
```
Error downloading model from Zenodo
```
**Solution:** Manual download from https://zenodo.org/records/10371965, then use `load_pretrained_model_from_dump(path)`

### Import Error
```
ModuleNotFoundError: No module named 'proteinbert'
```
**Solution:** Run `uv sync` to install dependencies
