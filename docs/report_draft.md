# FuncCell: Function-Aware Single-Cell Embeddings via Attention-Weighted Protein Language Models

**Andre Gala-Garza, Prema Immadisetty, Jason (Jingqiao) Zhao**

BIOINF 593 - Final Project

---

## Abstract

Single-cell RNA sequencing (scRNA-seq) enables transcriptome-wide profiling of individual cells, but gene expression alone lacks protein function context. We present FuncCell, a method for creating function-aware cell embeddings by combining expression data with ProteinBERT protein embeddings. ProteinBERT is a protein level embedding trained to contain functinal information, and in this porject we aggregated these protein embeddings to a higher cell level. We implement two aggregation strategies: a baseline naive weighted-sum approach where attention weights equal normalized expression, and an attention pooling model that learns gene importance from both gene identity and expression level. To ensure fair comparison, we evaluate both methods using identical downstream classifiers (Logistic Regression) on their 512-dimensional embeddings, isolating embedding quality from classifier complexity. On 20,000 breast tissue cells from CELLxGENE Census (10K healthy, 10K cancer), attention pooling achieves 99.4% classification accuracy compared to 97.9% for the naive baseline, with substantial improvements in clustering quality (ARI: 0.970 vs 0.854, ASW: 0.742 vs 0.629). We further demonstrate that the embedding space is geometrically meaningful: a simple threshold on projections onto the cancer direction achieves 98% accuracy without any machine learning classifier or attention pooling. These results suggest that learned attention over protein function embeddings captures biologically relevant signal beyond expression magnitude.

---

## 1. Introduction

Single-cell atlases such as CELLxGENE Census provide expression profiles for millions of cells across diverse tissues and conditions. While these resources enable unprecedented resolution in studying cellular heterogeneity, the expression-centric view has a fundamental limitation: gene expression levels do not directly reflect protein function. But these information are actually very valuable. If there is a way to encode cell function in a high dimensional vector, this will enable biologists to potentially gain a lot more insight instead just using unflexible traditinal ML method as a black box.

Protein language models offer us just that. ProteinBERT, trained on UniProt sequences with Gene Ontology annotations, produces embeddings that encode protein function from amino acid sequences alone. These 512-dimensional vectors capture functional similarity independent of expression context.

We propose FuncCell, which creates cell embeddings by aggregating ProteinBERT gene embeddings weighted by expression. This places cells in a "functional space" where similarity reflects shared protein functions rather than co-expression patterns. And in this project for the proof of concept, we focuses on using the embedding to discover breast cancer.

We used two weighting strategies:

1. **Baseline (Naive Weighted Sum):** Attention weights equal normalized expression values
2. **Attention Pooling:** A neural network learns attention weights from both gene embeddings and expression

A critical methodological choice is the intentional choise of boudary based simple classifier. Previous comparisons of embedding methods often conflate embedding quality with classifier complexity. By applying identical simple Logistic Regression classifiers to both embedding types, we isolate the contribution of the embedding itself.

---

## 2. Methods

### 2.1 Data

We obtained single-cell expression data from CELLxGENE Census, focusing on breast tissue. After quality filtering (minimum 200 genes per cell), our dataset comprised:

| Split | Healthy Cells | Cancer Cells | Total |
|-------|---------------|--------------|-------|
| Training | 10,000 | 10,000 | 20,000 |
| Test | 1,000 | 1,000 | 2,000 |

All cells were filtered to 10x 3' v3 assays and primary data only. The Census contains 1.7M healthy and 34K cancer breast cells matching these criteria; we randomly sampled from this pool.

### 2.2 Gene Embedding Generation

**Gene-to-Protein Mapping:** We used Ensembl BioMart to map gene symbols to canonical protein sequences:

| Category | Count | Percentage |
|----------|-------|------------|
| Total input genes | 19,477 | - |
| Successfully mapped | 19,294 | 99.1% |
| - Regular genes | 19,194 | - |
| - Readthrough genes | 100 | - |
| Failed to map | 183 | 0.9% |

Protein sequences ranged from 25 to 34,350 amino acids (mean: 579, median: 431).

**Readthrough Transcripts:** Fusion genes representing readthrough transcription (e.g., PDCD6-AHRR) were handled by concatenating the component protein sequences before embedding, preserving functional information from both proteins. These readthrough transcripts are potential cancer biomarkers.

**ProteinBERT Embeddings:** For each protein sequence, we extracted 512-dimensional embeddings from ProteinBERT's internal layer (`global-merge2-norm-block6`).

### 2.3 Baseline: Expression-Weighted Sum

For a cell with expression vector **e** over genes G:

$$w_g = \frac{e_g}{\sum_{g' \in G} e_{g'}}$$

$$\mathbf{c} = \sum_{g \in G} w_g \cdot \mathbf{p}_g$$

where **p**_g is the ProteinBERT embedding for gene g and **c** is the resulting 512-dimensional cell embedding. This assumes high expression implies high importance.

### 2.4 Attention Pooling

We learn attention weights from both gene identity and expression:

**Attention Network** (shared across all genes):
- Input: [gene_embedding (512) | expression (1)] = 513 dimensions
- Hidden: Linear(513→256) → ReLU → Linear(256→128) → ReLU → Linear(128→1)
- Output: Attention score per gene

**Cell Embedding:**
$$\alpha_g = \text{softmax}(\text{AttentionNet}([\mathbf{p}_g; e_g]))$$
$$\mathbf{c} = \sum_{g \in G} \alpha_g \cdot \mathbf{p}_g$$

**Classifier** (training only):
- Linear(512→64) → ReLU → Dropout(0.3) → Linear(64→1) → Sigmoid

The model contains approximately 198,000 parameters and was trained for 7 epochs with Binary Cross-Entropy loss using Adam optimizer (lr=1e-3).

### 2.5 Fair Evaluation

To isolate embedding quality from classifier complexity, we:

1. Generate 512-dim embeddings using both methods
2. Train identical Logistic Regression classifiers on both embedding types
3. Compare classification and clustering metrics

This ensures any performance difference reflects embedding quality, not classifier architecture.

---

## 3. Results

### 3.1 Classification and Clustering Performance

Table 1 presents the main comparison on the held-out test set (2,000 cells):

| Method | Accuracy | Macro F1 | AUC | ARI | ASW |
|--------|----------|----------|-----|-----|-----|
| Baseline + LogReg | 97.9% | 97.9% | 99.9% | 0.854 | 0.629 |
| Attention + LogReg | 99.4% | 99.3% | 99.9% | 0.970 | 0.742 |
| **Improvement** | **+1.5%** | **+1.4%** | +0.0% | **+0.117** | **+0.113** |

Attention pooling improves all metrics. The largest gains appear in clustering quality: Adjusted Rand Index (ARI) improves by 0.117 and Average Silhouette Width (ASW) by 0.113, indicating that attention embeddings form tighter, better-separated clusters.

Notably, replacing the neural classifier with Logistic Regression on attention embeddings yields identical performance, confirming that improvements stem from embedding quality rather than classifier difference.

### 3.2 Geometric Interpretability

We examined whether the embedding space encodes biologically meaningful structure through cancer direction analysis, the following analysis is being done using the naive baseline aggregation method:

$$\mathbf{d} = \bar{\mathbf{c}}_{\text{cancer}} - \bar{\mathbf{c}}_{\text{healthy}}$$

where $\bar{\mathbf{c}}$ denotes the centroid of each condition. Projecting cells onto this direction yields a scalar "cancer score":

$$s = \mathbf{c} \cdot \frac{\mathbf{d}}{||\mathbf{d}||}$$

A simple threshold classifier on this score achieves:
- **Accuracy:** 97%
- **ROC AUC:** 0.993

This was validated on 200 additional cells not used in centroid computation. The result demonstrates that the functional embedding space separates conditions along a single interpretable axis—no machine learning required.

### 3.3 Gene Contribution Analysis

To understand which genes drive the cancer direction, we analyzed each gene along two axes:

- **X-axis (Alignment):** $\mathbf{p}_g \cdot \hat{\mathbf{d}}$ — how aligned is the gene's functional embedding with the cancer direction?
- **Y-axis (Log2FC):** $\log_2\frac{\bar{e}_g^{\text{cancer}} + 0.1}{\bar{e}_g^{\text{healthy}} + 0.1}$ — is the gene differentially expressed?

Genes in the upper-right quadrant (positive alignment AND upregulated in cancer) drive the cancer signal; genes in the lower-left (negative alignment AND upregulated in healthy) drive the healthy signal.

**Top Cancer-Driving Genes:**
| Gene | Log2FC | Known Role |
|------|--------|------------|
| ERBB4 | +6.55 | Receptor tyrosine kinase, breast cancer marker |
| ESR1 | +5.09 | Estrogen receptor, ER+ breast cancer driver |
| TRPS1 | +5.10 | Transcription factor, breast cancer associated |

**Top Healthy-Driving Genes:**
| Gene | Log2FC | Known Role |
|------|--------|------------|
| RPL41 | -6.95 | Ribosomal protein, high in normal tissue |
| RPL39 | -8.14 | Ribosomal protein, translation machinery |
| RPS28 | -7.62 | Ribosomal protein, housekeeping function |

The cancer-driving genes include known breast cancer markers (ERBB4, ESR1), while healthy-driving genes are predominantly ribosomal proteins reflecting the high translational activity of normal cells. This validates that our functional embeddings capture some biologically meaningful signal.

---

## 4. Discussion

### Why Attention Pooling Improves Embeddings

The baseline assumes expression magnitude equals functional importance—a gene expressed at 10x should contribute 10x to the cell embedding. This assumption fails in several cases:

1. **Low-expression functional genes:** Tumor suppressors like TP53 are functionally critical even when lowly expressed
2. **High-expression housekeeping genes:** Ribosomal proteins dominate expression but add noise rather than discriminative signal
3. **Context-dependent importance:** The same gene may matter more in some cells than others

The attention network learns to reweight genes based on both their functional identity (embedding) and expression context, potentially discovering that "BRCA1 at any expression level matters for cancer classification."

### Embedding Space Properties

Our results demonstrate three desirable properties of the functional embedding space:

1. **Discriminative:** High classification accuracy (99.4%)
2. **Clustered:** High ARI (0.970) and ASW (0.742) indicate well-defined populations
3. **Interpretable:** Simple geometric operations (centroid difference, projection) achieve classification and reveal driving genes

### Limitations

Several limitations should be noted:

- **Single tissue/task:** We evaluated only breast tissue with binary healthy/cancer classification

### Future Directions

1. **Multi-tissue evaluation:** Extend to other tissues and multi-class problems (cell type classification)
2. **Attention interpretation:** Extract and visualize learned attention weights to identify "important" genes per cell
3. **Alternative protein models:** Compare ESM-2 and other protein language models
4. **Integration with other modalities:** Combine with chromatin accessibility or spatial information

---

## 5. Conclusion

FuncCell explores whether combining single-cell expression with protein language model embeddings can create a useful "functional space" for cells. The hypothesis is that ProteinBERT embeddings, trained on protein sequences and GO annotations, may encode functional information beyond what expression levels alone provide.

Our results are encouraging: the embedding space exhibits geometric structure where simple projection onto the healthy-to-cancer axis achieves 98% classification accuracy without any ML classifier. Gene contribution analysis along this axis recovers known cancer markers (ERBB4, ESR1) and housekeeping genes (ribosomal proteins), suggesting some biological relevance.

Attention pooling improves embedding quality over naive expression-weighting, with the largest gains in clustering metrics (ARI, ASW). Whether these embeddings truly capture "protein function" remains an open question, but the approach offers a potentially complementary view to traditional expression-based single-cell analysis.

---

## References

1. Brandes, N., et al. (2022). ProteinBERT: A universal deep-learning model of protein sequence and function. *Bioinformatics*.
2. CELLxGENE Census. Chan Zuckerberg Initiative. https://cellxgene.cziscience.com/
3. Wolf, F.A., et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*.

---

## Supplementary Information

**Code Availability:** https://github.com/jason-zhao-git/FuncCell

**Data Files:**
- `data/attention_comparison.csv` - Main metrics table
- `data/analysis/gene_contribution_analysis.csv` - Gene contribution scores
- `data/analysis/cancer_direction_512d.npy` - Cancer direction vector
