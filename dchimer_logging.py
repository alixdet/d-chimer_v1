#!/usr/bin/python3
# -*-coding:utf-8 -*

"""
Logging configuration for d-chimer.

Provides a centralized logging setup for all d-chimer modules.
"""

import logging
import sys
from pathlib import Path

def setup_logging(log_level=logging.INFO, log_file=None):
    """
    Configure logging for d-chimer.

    Args:
        log_level: logging level (default: logging.INFO)
        log_file: optional path to log file. If None, logs only to stdout.

    Returns:
        logger: configured logger instance
    """
    logger = logging.getLogger('dchimer')
    logger.setLevel(log_level)

    # Avoid adding handlers multiple times
    if logger.handlers:
        return logger

    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Add stdout handler
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(log_level)
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)

    # Add file handler if specified
    if log_file:
        try:
            log_path = Path(log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(log_level)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        except (OSError, IOError) as e:
            logger.warning(f"Could not create log file {log_file}: {e}")

    return logger


def get_logger(name=None):
    """
    Get a logger instance for a specific module.

    Args:
        name: module name

    Returns:
        logger: logger instance
    """
    if name is None:
        name = 'dchimer'
    return logging.getLogger(name)
