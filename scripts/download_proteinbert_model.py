#!/usr/bin/env python3
"""
Download ProteinBERT pretrained model from Zenodo.

This script downloads the model to the local project directory
to avoid the interactive FTP prompt and connection issues.

Usage:
    python scripts/download_proteinbert_model.py
"""

import sys
import urllib.request
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Model directory in the project
MODEL_DIR = project_root / "models" / "proteinbert"
MODEL_PATH = MODEL_DIR / "full_go_epoch_92400_sample_23500000.pkl"

# Zenodo download URL (more reliable than FTP)
ZENODO_URL = "https://zenodo.org/records/10371965/files/full_go_epoch_92400_sample_23500000.pkl?download=1"


def download_model():
    """Download ProteinBERT model from Zenodo."""

    # Check if model already exists
    if MODEL_PATH.exists():
        size_mb = MODEL_PATH.stat().st_size / (1024 * 1024)
        print(f"✓ Model already exists at {MODEL_PATH}")
        print(f"  Size: {size_mb:.1f} MB")
        return True

    # Create model directory
    MODEL_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Downloading ProteinBERT model from Zenodo...")
    print(f"  Source: {ZENODO_URL}")
    print(f"  Target: {MODEL_PATH}")
    print(f"  Size: ~500 MB (this may take a few minutes)")
    print()

    try:
        # Download with progress reporting
        def reporthook(block_num, block_size, total_size):
            downloaded = block_num * block_size
            percent = min(100, downloaded * 100 / total_size)
            mb_downloaded = downloaded / (1024 * 1024)
            mb_total = total_size / (1024 * 1024)
            print(f"\r  Progress: {percent:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)", end='')

        urllib.request.urlretrieve(ZENODO_URL, MODEL_PATH, reporthook)
        print()  # New line after progress

        # Verify download
        if MODEL_PATH.exists():
            size_mb = MODEL_PATH.stat().st_size / (1024 * 1024)
            print(f"\n✓ Model downloaded successfully!")
            print(f"  Location: {MODEL_PATH}")
            print(f"  Size: {size_mb:.1f} MB")
            return True
        else:
            print(f"\n✗ Download failed: File not found")
            return False

    except Exception as e:
        print(f"\n✗ Download failed: {e}")
        print("\nAlternative: Download manually from:")
        print(f"  {ZENODO_URL}")
        print(f"And save to: {MODEL_PATH}")
        return False


if __name__ == "__main__":
    success = download_model()
    sys.exit(0 if success else 1)
