"""File I/O utilities for funcCell project."""

import pickle
from pathlib import Path
from typing import Any
import logging

logger = logging.getLogger(__name__)


def ensure_dir(path: Path) -> Path:
    """
    Ensure directory exists, create if it doesn't.

    Args:
        path: Directory path

    Returns:
        Path object
    """
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_pickle(obj: Any, filepath: Path) -> None:
    """
    Save object to pickle file.

    Args:
        obj: Object to save
        filepath: Path to save file
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, 'wb') as f:
        pickle.dump(obj, f)
    logger.info(f"Saved pickle to: {filepath}")


def load_pickle(filepath: Path) -> Any:
    """
    Load object from pickle file.

    Args:
        filepath: Path to pickle file

    Returns:
        Loaded object
    """
    with open(filepath, 'rb') as f:
        obj = pickle.load(f)
    logger.info(f"Loaded pickle from: {filepath}")
    return obj


def save_text(text: str, filepath: Path) -> None:
    """
    Save text to file.

    Args:
        text: Text content
        filepath: Path to save file
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, 'w') as f:
        f.write(text)
    logger.info(f"Saved text to: {filepath}")


def save_list(items: list, filepath: Path) -> None:
    """
    Save list of items to text file (one per line).

    Args:
        items: List of items
        filepath: Path to save file
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, 'w') as f:
        for item in items:
            f.write(f"{item}\n")
    logger.info(f"Saved {len(items)} items to: {filepath}")


def load_list(filepath: Path) -> list:
    """
    Load list from text file (one per line).

    Args:
        filepath: Path to file

    Returns:
        List of items
    """
    with open(filepath, 'r') as f:
        items = [line.strip() for line in f if line.strip()]
    logger.info(f"Loaded {len(items)} items from: {filepath}")
    return items
