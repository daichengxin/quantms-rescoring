"""
Centralized logging configuration for quantmsrescore.

This module provides a consistent logging setup across the entire package,
with customizable log levels and formatters.
"""

import logging
import sys
import warnings
from typing import Optional


def configure_logging(log_level: str = "INFO") -> None:
    """
    Configure the logging system for the quantmsrescore package.

    This function sets up a consistent logging configuration with the specified
    log level and a standard format. It also suppresses all warnings from the
    pyopenms module.

    Parameters
    ----------
    log_level : str, optional
        The logging level to use (e.g., "DEBUG", "INFO", "WARNING", "ERROR").
        Default is "INFO".
    """
    # Convert string log level to numeric value
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level}")

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(numeric_level)

    # Clear any existing handlers to avoid duplicate logging
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Create console handler with a standard formatter
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(numeric_level)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    # Suppress all warnings from pyopenms
    warnings.filterwarnings("ignore", module="pyopenms")


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """
    Get a logger with the specified name.

    This function returns a logger with the specified name, which inherits
    the configuration from the root logger. If no name is provided, the root
    logger is returned.

    Parameters
    ----------
    name : str, optional
        The name of the logger to get. If None, the root logger is returned.

    Returns
    -------
    logging.Logger
        A logger with the specified name.
    """
    return logging.getLogger(name)