"""Logging utilities for funcCell project."""

import logging
import sys
from pathlib import Path
from datetime import datetime


def setup_logger(name: str, log_file: Path | None = None, level: int = logging.INFO) -> logging.Logger:
    """
    Setup logger with console and optional file handlers.

    Args:
        name: Logger name
        log_file: Optional path to log file
        level: Logging level (default: INFO)

    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Remove existing handlers
    logger.handlers = []

    # Console handler with formatting
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_format = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)

    # File handler if specified
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(console_format)
        logger.addHandler(file_handler)

    return logger


def log_step(logger: logging.Logger, step_name: str):
    """
    Decorator to log function execution time.

    Args:
        logger: Logger instance
        step_name: Name of the step being logged
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            logger.info(f"{'='*60}")
            logger.info(f"Starting: {step_name}")
            logger.info(f"{'='*60}")
            start_time = datetime.now()

            try:
                result = func(*args, **kwargs)
                elapsed = (datetime.now() - start_time).total_seconds()
                logger.info(f"✓ Completed: {step_name} in {elapsed:.2f}s")
                return result
            except Exception as e:
                elapsed = (datetime.now() - start_time).total_seconds()
                logger.error(f"✗ Failed: {step_name} after {elapsed:.2f}s")
                logger.error(f"Error: {str(e)}")
                raise

        return wrapper
    return decorator
