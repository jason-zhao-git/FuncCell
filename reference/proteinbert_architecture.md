# ProteinBERT Model Architecture

**Date:** 2025-11-21
**Model:** ProteinBERT (Brandes et al., Bioinformatics 2022)
**Source:** nadavbra/protein_bert
**Discovered through:** Direct model inspection and testing

---

## Executive Summary

ProteinBERT is a protein language model with a **dual-pathway architecture**: local (per-amino-acid) and global (whole-protein) representations. The model outputs are task-specific predictions (GO terms, amino acid types), NOT general-purpose embeddings. The true 512-dimensional protein embeddings must be extracted from an internal layer.

---

## Model Structure Overview

### Input Layer
- **input-seq**: (None, seq_len) - Amino acid sequence tokens
- **input-annotations**: (None, 8943) - GO annotation context (optional)

### Core Architecture
- **6 transformer-like blocks** with global attention mechanism
- Each block has two parallel pathways:
  1. **Local pathway**: Per-amino-acid representations (128-dim)
  2. **Global pathway**: Whole-protein representation (512-dim)

### Output Heads (Task-Specific)
- **output-seq**: (None, seq_len, 26) - Per-amino-acid predictions (20 amino acids + 6 special tokens)
- **output-annotations**: (None, 8943) - GO term predictions (multi-label)

---

## Layer-by-Layer Structure

### Block Pattern (Repeated 6 times)

Each block (numbered 1-6) contains:

```
Block N:
├── Global Pathway:
│   ├── global-to-seq-dense-blockN      (None, 128)  - Project global to sequence dim
│   ├── global-to-seq-reshape-blockN    (None, 1, 128) - Reshape for broadcasting
│   ├── global-dense1-blockN            (None, 512)  - Feed-forward layer 1
│   ├── global-attention-blockN         (None, 512)  - Aggregate from local pathway
│   ├── global-merge1-blockN            (None, 512)  - Residual connection
│   ├── global-merge1-norm-blockN       (None, 512)  - Layer normalization
│   ├── global-dense2-blockN            (None, 512)  - Feed-forward layer 2
│   ├── global-merge2-blockN            (None, 512)  - Residual connection
│   └── global-merge2-norm-blockN       (None, 512)  - Layer normalization ← **Extract here!**
│
└── Local Pathway:
    ├── narrow-conv-blockN              (None, seq_len, 128) - Narrow convolution
    ├── wide-conv-blockN                (None, seq_len, 128) - Wide convolution
    ├── seq-merge1-blockN               (None, seq_len, 128) - Merge convolutions
    ├── seq-merge1-norm-blockN          (None, seq_len, 128) - Layer normalization
    ├── seq-dense-blockN                (None, seq_len, 128) - Feed-forward
    ├── seq-merge2-blockN               (None, seq_len, 128) - Residual connection
    └── seq-merge2-norm-blockN          (None, seq_len, 128) - Layer normalization
```

### Key Layers for Embedding Extraction

**IMPORTANT:** The 512-dim global representation is at:
```
global-merge2-norm-block6  (None, 512)  ← Final global representation
```

This is the **true protein embedding** layer, NOT the output heads!

---

## Complete Layer List (102 layers total)

### Input & Embedding Layers (0-4)
```
0.  input-annotations                   (None, 8943)
1.  input-seq                           (None, seq_len)
2.  dense-global-input                  (None, 512)         - Initial GO embedding
3.  embedding-seq-input                 (None, seq_len, 128) - Amino acid embeddings
```

### Block 1 (5-19)
```
4.  global-to-seq-dense-block1          (None, 128)
5.  global-to-seq-reshape-block1        (None, 1, 128)
6.  narrow-conv-block1                  (None, seq_len, 128)
7.  wide-conv-block1                    (None, seq_len, 128)
8.  seq-merge1-block1                   (None, seq_len, 128)
9.  seq-merge1-norm-block1              (None, seq_len, 128)
10. seq-dense-block1                    (None, seq_len, 128)
11. seq-merge2-block1                   (None, seq_len, 128)
12. seq-merge2-norm-block1              (None, seq_len, 128)
13. global-dense1-block1                (None, 512)
14. global-attention-block1             (None, 512)
15. global-merge1-block1                (None, 512)
16. global-merge1-norm-block1           (None, 512)
17. global-dense2-block1                (None, 512)
18. global-merge2-block1                (None, 512)
19. global-merge2-norm-block1           (None, 512)
```

### Blocks 2-5 (20-83)
Similar structure repeated for blocks 2, 3, 4, 5...

### Block 6 (84-99) - FINAL BLOCK
```
84. global-to-seq-dense-block6          (None, 128)
85. global-to-seq-reshape-block6        (None, 1, 128)
86. narrow-conv-block6                  (None, seq_len, 128)
87. wide-conv-block6                    (None, seq_len, 128)
88. seq-merge1-block6                   (None, seq_len, 128)
89. seq-merge1-norm-block6              (None, seq_len, 128)
90. seq-dense-block6                    (None, seq_len, 128)
91. seq-merge2-block6                   (None, seq_len, 128)
92. seq-merge2-norm-block6              (None, seq_len, 128)
93. global-dense1-block6                (None, 512)
94. global-attention-block6             (None, 512)
95. global-merge1-block6                (None, 512)
96. global-merge1-norm-block6           (None, 512)
97. global-dense2-block6                (None, 512)
98. global-merge2-block6                (None, 512)
99. global-merge2-norm-block6           (None, 512)  ← **EXTRACT THIS FOR EMBEDDINGS!**
```

### Output Heads (100-101)
```
100. output-seq                         (None, seq_len, 26)    - Amino acid predictions
101. output-annotations                 (None, 8943)           - GO term predictions
```

---

## Dimension Summary

| Representation | Dimension | Location | Purpose |
|----------------|-----------|----------|---------|
| **Local (per-AA)** | 128 | `seq-merge2-norm-block*` | Per-amino-acid features |
| **Global (protein)** | 512 | `global-merge2-norm-block6` | **Whole-protein embedding** ✓ |
| AA Predictions | 26 | `output-seq` | Task-specific (NOT embedding) |
| GO Predictions | 8943 | `output-annotations` | Task-specific (NOT embedding) |

---

## Common Mistakes & Corrections

### ❌ WRONG: Using Model Outputs as Embeddings

```python
# These are PREDICTIONS, not embeddings!
model = pretrained_model_generator.create_model(seq_len)
outputs = model.predict(encoded_x)

# WRONG: Using GO predictions
embedding = outputs[1]  # Shape: (1, 8943) ← GO predictions, NOT embeddings!
```

### ❌ WRONG: Using Concatenated Hidden Layers

```python
from proteinbert.conv_and_global_attention_model import get_model_with_hidden_layers_as_outputs

# This concatenates ALL hidden layers + GO predictions
model = get_model_with_hidden_layers_as_outputs(base_model)
local, global_concat = model.predict(encoded_x)

# WRONG: This is 15,599-dim (8943 GO + 6656 hidden states)
embedding = global_concat[0]  # Shape: (15599,) ← Concatenated outputs, NOT clean embeddings!

# Where 15,599 comes from:
# - 8,943 (GO term predictions)
# - 6,656 (13 layers × 512 dim per layer)
# Total: 15,599 dimensions
```

### ✓ CORRECT: Extracting Internal Global Layer

```python
import tensorflow as tf
from proteinbert import load_pretrained_model

# Load model
pretrained_model_generator, input_encoder = load_pretrained_model()
base_model = pretrained_model_generator.create_model(seq_len=2048)

# Extract the 512-dim global representation layer
final_global_layer = base_model.get_layer('global-merge2-norm-block6')

# Create new model that outputs the embedding
embedding_model = tf.keras.Model(
    inputs=base_model.inputs,
    outputs=final_global_layer.output
)

# Get clean 512-dim embedding
encoded_x = input_encoder.encode_X([sequence], seq_len)
embedding = embedding_model.predict(encoded_x)[0]  # Shape: (512,) ✓ CORRECT!
```

---

## Model Input Format

### Sequence Encoding
```python
encoded_x = input_encoder.encode_X(sequences, seq_len)
# Returns tuple: (seq_array, annotation_array)
# - seq_array: (batch, seq_len) - Tokenized amino acids
# - annotation_array: (batch, 8943) - GO context (zeros if not provided)
```

### Token Vocabulary
- 20 standard amino acids: ACDEFGHIKLMNPQRSTVWY
- Special tokens: START, END, PAD, MASK, etc.
- Total: 26 tokens

---

## Global Attention Mechanism

ProteinBERT uses a **linear complexity attention** (O(n) instead of O(n²)):

1. **Local representations** attend within neighborhoods
2. **Global representation** aggregates from all positions via attention
3. **Global-to-local broadcast** injects global context back into sequence
4. This enables processing of very long sequences (up to 34,350 amino acids)

### Key Innovation
Standard transformer attention:
```
Complexity: O(n²) where n = sequence length
Max sequence: ~1000 amino acids
```

ProteinBERT global attention:
```
Complexity: O(n) where n = sequence length
Max sequence: 34,000+ amino acids ✓
```

---

## Pretraining Details

### Dataset
- **Source:** UniRef90
- **Size:** 106 million protein sequences
- **Training time:** 28 days on single GPU
- **Throughput:** 280 proteins/second during training

### Pretraining Tasks
1. **Masked Language Modeling:** Predict masked amino acids
2. **GO Term Prediction:** Multi-label classification for Gene Ontology terms

### Why GO Terms Matter
The model was pretrained to predict GO terms, which forces it to learn functional representations. The GO prediction head (8943-dim) is task-specific and should NOT be used as a general-purpose embedding. Instead, extract the 512-dim global representation that the model uses internally.

---

## Usage Recommendations

### For General Protein Embeddings
✓ Extract `global-merge2-norm-block6` (512-dim)
- Clean, general-purpose protein representation
- Learned functional properties from 106M proteins
- Suitable for downstream tasks

### For GO Term Prediction
✓ Use `output-annotations` (8943-dim)
- Task-specific predictions
- Multi-label probabilities for GO terms
- Good for functional annotation

### For Per-Amino-Acid Analysis
✓ Use `output-seq` (seq_len × 26)
- Per-residue predictions
- Amino acid type probabilities
- Good for sequence analysis

---

## Implementation in funcCell

```python
# File: src/model/proteinbert_embeddings.py

class ProteinBERTEmbedder:
    def __init__(self, seq_len: int = 2048):
        # Load base model
        self.pretrained_model_generator, self.input_encoder = load_pretrained_model(...)
        base_model = self.pretrained_model_generator.create_model(seq_len)

        # Extract 512-dim global layer
        final_global_layer = base_model.get_layer('global-merge2-norm-block6')

        # Create embedding model
        self.model = tf.keras.Model(
            inputs=base_model.inputs,
            outputs=final_global_layer.output  # (None, 512)
        )

    def embed_sequence(self, sequence: str) -> np.ndarray:
        encoded_x = self.input_encoder.encode_X([sequence], self.seq_len)
        embedding = self.model.predict(encoded_x, batch_size=1, verbose=0)
        return embedding[0]  # Shape: (512,)
```

---

## Verification

### Test Code
```python
# Verify embedding dimension
embedder = ProteinBERTEmbedder(seq_len=2048)
test_seq = "MKTIIALSYIFCLVFA"
embedding = embedder.embed_sequence(test_seq)

assert embedding.shape == (512,), f"Expected (512,) but got {embedding.shape}"
assert not np.isnan(embedding).any(), "Embedding contains NaN values"
assert not np.isinf(embedding).any(), "Embedding contains Inf values"

print(f"✓ Embedding shape: {embedding.shape}")
print(f"✓ Embedding magnitude: {np.linalg.norm(embedding):.2f}")
print(f"✓ First 10 values: {embedding[:10]}")
```

### Expected Output
```
✓ Embedding shape: (512,)
✓ Embedding magnitude: 8.45
✓ First 10 values: [-0.02864169 -0.12259582 -0.05216378 -0.07679206  0.23831926
 -0.07313731 -0.05246332 -0.36069673 -0.14692229 -0.20581266]
```

---

## References

- **Paper:** Brandes, N., Ofer, D., Peleg, Y., Rappoport, N., & Linial, M. (2022). ProteinBERT: a universal deep-learning model of protein sequence and function. *Bioinformatics*, 38(8), 2102-2110.
- **GitHub:** https://github.com/nadavbra/protein_bert
- **Model Weights:** https://zenodo.org/records/10371965
- **Hugging Face:** https://huggingface.co/GrimSqueaker/proteinBERT

---

## Discovery Notes

**Date:** 2025-11-21

Initial confusion stemmed from multiple sources claiming different embedding dimensions:
- Generic BERT documentation: 1024-dim
- Model output shape: 8943-dim (GO predictions)
- Hidden layer concatenation: 15,599-dim (all layers + GO)

**Resolution:** Direct model inspection revealed the true 512-dim global representation at `global-merge2-norm-block6`, which is the clean protein embedding learned during pretraining, separate from task-specific prediction heads.

This layer contains the learned functional representation that generalizes across proteins, making it suitable for downstream tasks like cancer cell classification in funcCell.
