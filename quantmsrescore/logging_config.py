"""
Centralized logging configuration for quantmsrescore.

This module provides a consistent logging setup across the entire package,
with customizable log levels and formatters.
"""

import logging
import sys
import warnings
import io
import re
from typing import Optional

class IgnoreSpecificWarnings(logging.Filter):
    def filter(self, record):
        # Check for any isotope-related atom warnings
        message = record.getMessage()
        if "Could not add the following atom:" in message:
            return False
        return True

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
    warnings.filterwarnings("ignore", module="ms2pip")
    warnings.filterwarnings("ignore", module="ms2rescore")
    warnings.filterwarnings("ignore", module="xgboost")
    warnings.filterwarnings("ignore", module="tensorflow")
    warnings.filterwarnings("ignore", module="deeplc")

    # Ignore annoying warning from ms2pip
    root_logger.addFilter(IgnoreSpecificWarnings())

    # Apply filter to all loggers, not just root
    for logger_name in logging.root.manager.loggerDict:
        logger = logging.getLogger(logger_name)
        logger.addFilter(IgnoreSpecificWarnings())
    
    # Suppress specific warnings using multiple approaches
    warnings.filterwarnings("ignore", message=".*Could not add the following atom.*")
    warnings.filterwarnings("ignore", message=".*\\[[0-9]+\\].*")  # Match any isotope notation like [13], [15], etc.

    # Reduce the log level for this specific warning pattern
    logging.getLogger("ms2pip").setLevel(logging.ERROR)

    # Capture warnings and redirect them to the logging system
    # This helps catch warnings that might bypass the regular filters
    original_showwarning = warnings.showwarning

    def custom_showwarning(message, category, filename, lineno, file=None, line=None):
        # Check if this is the specific warning we want to ignore
        msg_str = str(message)
        # Match any "Could not add the following atom" warning or any isotope notation
        if "Could not add the following atom" in msg_str or re.search(r'\[[0-9]+\]', msg_str):
            return  # Completely suppress the warning
        # For all other warnings, use the original handler
        return original_showwarning(message, category, filename, lineno, file, line)

    warnings.showwarning = custom_showwarning



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