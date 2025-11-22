# ProteinBERT Embedding Generation

**Date:** 2025-11-21 (Updated)
**Model:** ProteinBERT (pretrained on 106M proteins)
**Total Proteins:** 19,294 (19,194 regular + 100 readthrough)
**Embedding Dimension:** 512 (CORRECTED - was incorrectly 8943 or 15599)

## Summary

ProteinBERT embeddings provide protein function representations for all mapped genes. Each protein sequence (variable length: 25-34,350 amino acids) is compressed into a fixed **512-dimensional global representation** capturing functional properties learned from 106 million proteins.

## IMPORTANT: Correct Embedding Extraction (Updated 2025-11-21)

**Critical Discovery:** The ProteinBERT model outputs are NOT the protein embeddings. The model has multiple outputs:
- `output-seq`: (batch, seq_len, 26) - per-amino-acid predictions
- `output-annotations`: (batch, 8943) - GO term predictions (task-specific)

**The true 512-dim global protein embedding** is hidden inside the model at layer `global-merge2-norm-block6`.

### Incorrect Methods (DO NOT USE)
```python
# ❌ WRONG: Using GO predictions as embeddings
outputs = model.predict(encoded_x)
embedding = outputs[1]  # Shape: (1, 8943) - these are GO predictions, NOT embeddings!

# ❌ WRONG: Using concatenated hidden layers
from proteinbert.conv_and_global_attention_model import get_model_with_hidden_layers_as_outputs
model = get_model_with_hidden_layers_as_outputs(base_model)
local, global_concat = model.predict(encoded_x)
embedding = global_concat[0]  # Shape: (15599,) - concatenated GO + hidden states, NOT embeddings!
```

### Correct Method (USE THIS)
```python
import tensorflow as tf
from proteinbert import load_pretrained_model

# Load model
pretrained_model_generator, input_encoder = load_pretrained_model()
base_model = pretrained_model_generator.create_model(seq_len=2048)

# Extract the 512-dim global representation from internal layer
final_global_layer = base_model.get_layer('global-merge2-norm-block6')

# Create new model that outputs only the 512-dim embedding
embedding_model = tf.keras.Model(
    inputs=base_model.inputs,
    outputs=final_global_layer.output
)

# Get embeddings
encoded_x = input_encoder.encode_X([sequence], seq_len)
embedding = embedding_model.predict(encoded_x)[0]  # Shape: (512,) ✓ CORRECT!
```

## Model Architecture

### ProteinBERT Overview
- **Framework:** TensorFlow/Keras (not PyTorch)
- **Pretraining:** 106M proteins from UniRef90
- **Training time:** 28 days on GPU
- **Key innovation:** Global attention (linear complexity O(n) vs standard O(n²))
- **Advantage:** Can handle extremely long sequences up to 34k amino acids

### Dual Representation System
- **Local representations:** Per-amino acid embeddings (128-dim × sequence_length)
- **Global representations:** Whole-protein embeddings (512-dim per protein) ← **We use these**
- **Model outputs:** Task-specific predictions (26-dim per amino acid, 8943-dim GO terms) ← **NOT embeddings!**

## Implementation Details

### Configuration
```python
seq_len = 2048               # Maximum sequence length for encoding (UPDATED from 1024)
batch_size = 32              # Sequences processed per batch
readthrough_strategy = 'concat'  # Concatenate sequences (biologically accurate)
```

### Sequence Length Handling
- **Training lengths:** 128, 512, or 1024 tokens
- **Our choice:** 2048 tokens (UPDATED from 1024)
  - Covers ~95% of proteins (median=431aa, mean=579aa)
  - With 2048 buffer, even fewer proteins need truncation
  - Long sequences (>2048aa) are randomly subsampled
  - Model's global attention captures full-protein context from substrings

#### Truncation Strategy for Long Proteins (>2048aa)

**Dataset Statistics:**
- With 2048aa buffer, fewer proteins need truncation compared to 1024aa
- Longest: TTN (titin) at 34,350aa
- Other examples: MUC16 (14,507aa), SYNE1 (8,797aa), NEB (8,525aa)

**Chosen Approach: Single Random Substring**
```python
# For sequences >2048aa, randomly select one 2048aa window
if len(sequence) > 2048:
    start = random.randint(0, len(sequence) - 2048)
    sequence = sequence[start:start + 2048]
# Embed once → one 512-dim vector
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

**Solution:** Concatenate sequences then embed (UPDATED from mean pooling)
```python
# For readthrough gene with components (seq1, seq2)
# Concatenate into single fusion protein sequence
fusion_sequence = seq1 + seq2

# Embed the fusion as one protein
readthrough_embedding = embed_protein(fusion_sequence)  # 512-dim
```

**Why concatenation:**
- **Biologically accurate:** Readthrough transcripts produce continuous polypeptide chains
- Preserves dimensionality (all genes have 512-dim embeddings)
- Captures interactions between component proteins
- No information loss from averaging

**Alternative (no longer used): Mean pooling**
```python
# Old approach: embed separately then average
embedding1 = embed_protein(seq1)  # 512-dim
embedding2 = embed_protein(seq2)  # 512-dim
readthrough_embedding = (embedding1 + embedding2) / 2  # 512-dim
```
- This was simpler but loses biological accuracy
- Concatenation better reflects actual protein structure

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

### Model Loading (CORRECTED)
```python
import tensorflow as tf
from proteinbert import load_pretrained_model

# Load pretrained model (downloads ~500 MB on first run)
pretrained_model_generator, input_encoder = load_pretrained_model()

# Create base model
base_model = pretrained_model_generator.create_model(seq_len=2048)

# Extract 512-dim global representation from internal layer
final_global_layer = base_model.get_layer('global-merge2-norm-block6')

# Create embedding model
embedding_model = tf.keras.Model(
    inputs=base_model.inputs,
    outputs=final_global_layer.output
)
```

### Batch Processing
```python
# Encode sequences
encoded_x = input_encoder.encode_X(sequences, seq_len=2048)

# Get 512-dim global embeddings
global_embeddings = embedding_model.predict(encoded_x, batch_size=32)

# global_embeddings shape: (n_sequences, 512)
```

## Output

### File: `data/embeddings/gene_to_embedding.pkl`
```python
{
    'TP53': array([0.12, -0.34, ..., 0.56]),      # 512 floats
    'BRCA1': array([0.45, -0.12, ..., 0.89]),    # 512 floats
    'PDCD6-AHRR': array([0.23, 0.15, ..., 0.34]), # 512 floats (concatenated sequence)
    # ... 19,294 total genes
}
```

### Statistics
- **File size:** ~40 MB (19,294 × 512 × 4 bytes) - REDUCED from 660 MB!
- **Coverage:** 100% of mapped genes
- **Dimension:** All embeddings are 512-dim (CORRECTED from 8943 or 15599)
- **Quality:** No NaN or Inf values, all validated

## Validation

Post-generation checks performed:
1. **Coverage:** All 19,294 genes have embeddings
2. **Dimensionality:** All vectors are exactly 512-dim
3. **Validity:** No NaN or Inf values
4. **Magnitude:** Reasonable embedding norms (mean ± std logged)
5. **Readthrough:** All 100 fusion genes handled correctly (via concatenation)

## Usage Example

### Load Embeddings
```python
import pickle

with open('data/embeddings/gene_to_embedding.pkl', 'rb') as f:
    gene_to_embedding = pickle.load(f)

# Access embedding
tp53_embedding = gene_to_embedding['TP53']  # shape: (512,)
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
        cell_embeddings: Array of shape (n_cells, 512)
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

    # Get protein embeddings: (n_genes, 512)
    protein_embeddings = np.array([gene_to_embedding[g] for g in common_genes])

    # Weighted sum: (n_cells, 512)
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
- **Internal Representations:**
  - Local: (seq_len, 128) per-residue embeddings
  - **Global: (512,) whole-protein embedding** ← This is what we extract
- **Output Heads (task-specific, NOT embeddings):**
  - Sequence output: (seq_len, 26) per-amino-acid predictions
  - Annotation output: (8943,) GO term predictions

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
