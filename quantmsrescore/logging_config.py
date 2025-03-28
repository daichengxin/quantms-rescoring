"""
Centralized logging configuration for quantmsrescore.

This module provides a consistent logging setup across the entire package,
with customizable log levels and formatters.
"""

import logging
import sys
import warnings
import re
from typing import Optional


class IgnoreSpecificWarnings(logging.Filter):
    def filter(self, record):
        # Check for any warnings we want to ignore
        message = record.getMessage()

        # Isotope-related atom warnings
        if "Could not add the following atom:" in message:
            return False

        # DeepLC-related warnings
        if "Could not add the following value:" in message:
            return False

        if "Skipping the following (not in library):" in message:
            return False

        if "DeepLC tried to set intra op threads" in message:
            return False

        # OpenMS environment variable warning
        if re.search(r"Warning: OPENMS_DATA_PATH is not set", message):
            return False

        # CUDA and TensorFlow warnings
        if any(
            pattern in message
            for pattern in [
                "Unable to register cuDNN factory",
                "Unable to register cuBLAS factory",
                "computation placer already registered",
                "failed call to cuInit",
                "CUDA error",
            ]
        ):
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
    warnings.filterwarnings(
        "ignore",
        message="OPENMS_DATA_PATH environment variable already exists",
        category=UserWarning,
    )
    warnings.filterwarnings(
        action="ignore",
        message=".*OPENMS_DATA_PATH.*",
        category=UserWarning,
    )
    warnings.filterwarnings(
        "ignore",
        message="Warning: OPENMS_DATA_PATH environment variable already exists.*",
    )
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
    warnings.filterwarnings("ignore", message=".*Could not add the following value.*")
    warnings.filterwarnings("ignore", message=".*Skipping the following \(not in library\).*")
    warnings.filterwarnings("ignore", message=".*DeepLC tried to set intra op threads.*")
    warnings.filterwarnings(
        "ignore", message=".*\\[[0-9]+\\].*"
    )  # Match any isotope notation like [13], [15], etc.
    warnings.filterwarnings("ignore", message=".*OPENMS_DATA_PATH.*")

    # Suppress CUDA and TensorFlow warnings
    warnings.filterwarnings("ignore", message=".*Unable to register cuDNN factory.*")
    warnings.filterwarnings("ignore", message=".*Unable to register cuBLAS factory.*")
    warnings.filterwarnings("ignore", message=".*computation placer already registered.*")
    warnings.filterwarnings("ignore", message=".*failed call to cuInit.*")
    warnings.filterwarnings("ignore", message=".*CUDA error.*")
    # Reduce the log level for this specific warning pattern
    logging.getLogger("ms2pip").setLevel(logging.ERROR)

    # Capture warnings and redirect them to the logging system
    # This helps catch warnings that might bypass the regular filters
    original_showwarning = warnings.showwarning

    def custom_showwarning(message, category, filename, lineno, file=None, line=None):
        # Check if this is the specific warning we want to ignore
        msg_str = str(message)
        # Match any warnings we want to suppress
        cuda_tf_patterns = [
            "Unable to register cuDNN factory",
            "Unable to register cuBLAS factory",
            "computation placer already registered",
            "failed call to cuInit",
            "CUDA error",
        ]

        if (
            "Could not add the following atom" in msg_str
            or "Could not add the following value" in msg_str
            or "Skipping the following (not in library)" in msg_str
            or "DeepLC tried to set intra op threads" in msg_str
            or re.search(r"\[[0-9]+\]", msg_str)
            or "OPENMS_DATA_PATH" in msg_str
            or any(pattern in msg_str for pattern in cuda_tf_patterns)
        ):
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
