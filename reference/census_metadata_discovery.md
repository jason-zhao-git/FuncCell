# CELLxGENE Census Metadata Discovery Results

**Date:** 2025-11-18
**Census Version:** Latest
**Query:** Breast tissue data for cancer classification

## Summary

Discovered the correct metadata filter values for querying breast cancer vs. healthy tissue from CELLxGENE Census.

## Key Findings

### Correct Filter Fields

- **Tissue field:** Use `tissue_general == 'breast'` (NOT `tissue == 'breast tissue'`)
- **Disease field:** Use `disease` with exact values below
- **Assay field:** Use `assay == "10x 3' v3"` for standardization
- **Primary data:** Use `is_primary_data == True` to avoid duplicates

### Available Disease Values in Breast Tissue

| Disease | Cell Count |
|---------|-----------|
| `'normal'` | 5,433,737 |
| `'breast cancer'` | 641,674 |
| `'invasive ductal breast carcinoma'` | 601,153 |
| `'triple-negative breast carcinoma'` | 184,413 |
| `'estrogen-receptor positive breast cancer'` | 163,682 |
| `'invasive lobular breast carcinoma'` | 55,894 |
| `'breast carcinoma'` | 29,490 |
| `'HER2 positive breast carcinoma'` | 25,701 |
| `'invasive tubular breast carcinoma \|\| invasive lobular breast carcinoma'` | 17,024 |
| `'breast mucinous carcinoma'` | 11,370 |
| `'breast apocrine carcinoma'` | 8,232 |
| `'metaplastic breast carcinoma'` | 3,274 |

**Total breast tissue cells:** 7,175,644

### Available Assay Types in Breast Tissue

| Assay | Cell Count |
|-------|-----------|
| `'10x 3' v3'` | 5,091,015 ⭐ |
| `'10x 3' v2'` | 1,231,666 |
| `'10x 5' transcription profiling'` | 373,796 |
| `'10x 5' v1'` | 271,562 |
| `'10x 5' v2'` | 109,300 |
| `'10x multiome'` | 51,367 |
| `'inDrop'` | 38,219 |
| `'10x 3' transcription profiling'` | 7,709 |
| `'Smart-seq2'` | 1,010 |

### Tissue Values

| Tissue | Cell Count |
|--------|-----------|
| `'breast'` | 7,098,851 |
| `'upper outer quadrant of breast'` | 76,793 |

## Recommended Filters for funcCell

### For Healthy Breast Tissue:
```python
obs_value_filter = "tissue_general == 'breast' and disease == 'normal' and assay == \"10x 3' v3\" and is_primary_data == True"
```

**Expected cells:** ~1,726,582 cells

### For Breast Cancer Tissue:
```python
obs_value_filter = "tissue_general == 'breast' and disease == 'breast cancer' and assay == \"10x 3' v3\" and is_primary_data == True"
```

**Expected cells:** ~221,168 cells

### Alternative: Include Specific Carcinoma Types
```python
obs_value_filter = "tissue_general == 'breast' and disease in ['breast cancer', 'invasive ductal breast carcinoma', 'triple-negative breast carcinoma'] and assay == \"10x 3' v3\" and is_primary_data == True"
```

**Expected cells (carcinoma variants):** ~130,449 cells

## Benefits of These Filters

1. **Assay Standardization (`10x 3' v3`):**
   - 5M+ cells available (70% of all breast data)
   - Reduces technical batch effects
   - Ensures comparable gene detection across samples
   - Prevents model from learning assay-specific artifacts

2. **Primary Data Only (`is_primary_data == True`):**
   - Avoids duplicate/reprocessed cells
   - Ensures original count matrices
   - Cleaner training data for ML models

3. **Tissue General (`tissue_general`):**
   - Broader categorization than specific tissue
   - More consistent across datasets
   - Better for cross-study queries

## Data Quality Check

With these filters applied:
- ✅ **Healthy cells:** 1.7M cells (excellent for training)
- ✅ **Cancer cells:** 221K cells (sufficient for balanced dataset)
- ✅ **Ratio:** ~8:1 healthy:cancer (can downsample healthy for balance)
- ✅ **Technical consistency:** Single assay type (10x 3' v3)
- ✅ **No duplicates:** Primary data only

## Important: Census Schema Limitation

**Census does NOT have `feature_biotype` in var metadata!**
- Census only includes features with `feature_biotype='gene'`
- Cannot filter for protein-coding genes directly in Census query
- **Solution:** Use BioMart annotations (`data/raw/mart_export.txt`) to filter after querying

## Discovery Script

Run to rediscover metadata values:
```bash
uv run python scripts/discover_metadata.py
```

## References

- CELLxGENE Census Documentation: https://chanzuckerberg.github.io/cellxgene-census/
- Metadata Schema: https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md
- Collections:
  - Healthy: https://cellxgene.cziscience.com/collections/48259aa8-f168-4bf5-b797-af8e88da6637
  - Tumor: https://cellxgene.cziscience.com/collections/dea97145-f712-431c-a223-6b5f565f362a
