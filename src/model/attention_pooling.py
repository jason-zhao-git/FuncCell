"""Attention Pooling model for cell classification.

This module provides:
- AttentionPooling: Neural network that learns gene importance from embedding + expression
- CellDataset: PyTorch Dataset for loading cell data with gene embeddings
"""

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset
from scipy.sparse import issparse
import anndata


class AttentionPooling(nn.Module):
    """Attention-based pooling that learns gene importance from embedding + expression.

    Unlike simple weighted sum (where attention = expression), this model learns
    attention = f(embedding, expression), allowing it to discover which genes
    matter regardless of expression level.

    Args:
        embed_dim: Dimension of gene embeddings (default: 512 for ProteinBERT)
        hidden_dim: Hidden dimension of attention network (default: 128)
    """

    def __init__(self, embed_dim: int = 512, hidden_dim: int = 128):
        super().__init__()
        # Attention network: sees gene embedding (512) + expression (1) = 513
        self.attention = nn.Sequential(
            nn.Linear(embed_dim + 1, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, 1)
        )
        # Classifier head
        self.classifier = nn.Sequential(
            nn.Linear(embed_dim, 64),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )

    def forward(self, gene_embeddings, expression, expression_mask):
        """Forward pass.

        Args:
            gene_embeddings: (batch, n_genes, embed_dim) - frozen ProteinBERT embeddings
            expression: (batch, n_genes) - normalized expression values
            expression_mask: (batch, n_genes) - boolean mask for expressed genes

        Returns:
            pred: (batch, 1) - probability of cancer
            cell_embedding: (batch, embed_dim) - learned cell embedding
            attn_weights: (batch, n_genes) - attention weights per gene
        """
        # Concatenate expression to embeddings: (batch, n_genes, 513)
        gene_features = torch.cat([gene_embeddings, expression.unsqueeze(-1)], dim=-1)

        # Compute attention weights
        attn_scores = self.attention(gene_features).squeeze(-1)  # (batch, n_genes)
        attn_scores = attn_scores.masked_fill(~expression_mask, -1e9)
        attn_weights = torch.softmax(attn_scores, dim=-1)  # (batch, n_genes)

        # Weighted pool (still embed_dim output)
        cell_embedding = torch.bmm(attn_weights.unsqueeze(1), gene_embeddings).squeeze(1)

        # Classify
        return self.classifier(cell_embedding), cell_embedding, attn_weights


class CellDataset(Dataset):
    """Dataset for attention pooling model.

    Loads expression data from AnnData objects and pairs with gene embeddings.

    Args:
        adata_list: List of AnnData objects to concatenate
        labels_list: List of label arrays corresponding to each AnnData
        gene_to_embedding: Dict mapping gene names to embedding vectors
    """

    def __init__(self, adata_list, labels_list, gene_to_embedding):
        # Concatenate all adata
        self.adata = anndata.concat(adata_list, join='inner')
        self.labels = np.concatenate(labels_list)

        # Get common genes with embeddings
        self.gene_list = [g for g in self.adata.var_names if g in gene_to_embedding]
        print(f"Using {len(self.gene_list)} genes with embeddings")

        # Pre-compute gene embeddings matrix (same for all cells)
        self.gene_embeddings = np.stack(
            [gene_to_embedding[g] for g in self.gene_list]
        ).astype(np.float32)

        # Filter adata to common genes
        self.adata = self.adata[:, self.gene_list]

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        expr = self.adata[idx].X
        if issparse(expr):
            expr = expr.toarray().flatten()
        else:
            expr = np.asarray(expr).flatten()

        # Normalize expression (sum to 1)
        expr_norm = expr / (expr.sum() + 1e-10)

        mask = expr > 0
        return {
            'embeddings': torch.tensor(self.gene_embeddings, dtype=torch.float32),
            'expression': torch.tensor(expr_norm, dtype=torch.float32),
            'mask': torch.tensor(mask, dtype=torch.bool),
            'label': torch.tensor(self.labels[idx], dtype=torch.float32)
        }


def compute_baseline_embeddings(dataset):
    """Compute weighted sum embeddings (baseline method).

    This is the simple approach where attention = expression.

    Args:
        dataset: CellDataset instance

    Returns:
        labels: (n_cells,) array of labels
        embeddings: (n_cells, embed_dim) array of cell embeddings
    """
    embeddings = []
    labels = []

    for i in range(len(dataset)):
        item = dataset[i]
        expr = item['expression'].numpy()
        gene_emb = item['embeddings'].numpy()

        # Weighted sum: expression-weighted average of gene embeddings
        cell_emb = expr @ gene_emb  # (n_genes,) @ (n_genes, 512) = (512,)
        embeddings.append(cell_emb)
        labels.append(item['label'].item())

    return np.array(labels), np.stack(embeddings)
