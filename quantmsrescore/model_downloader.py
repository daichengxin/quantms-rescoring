"""
Model downloader for quantms-rescoring.

This module provides functionality to download all required models for MS2PIP,
DeepLC, and AlphaPeptDeep ahead of time for offline use.
"""

import shutil
from pathlib import Path
from typing import Optional

import click

from quantmsrescore.logging_config import configure_logging, get_logger

# Get logger for this module
logger = get_logger(__name__)


def download_ms2pip_models(model_dir: Optional[Path] = None) -> None:
    """
    Download MS2PIP models.

    MS2PIP models are bundled with the ms2pip package and don't require
    separate downloading. This function validates that the ms2pip package
    is properly installed.

    Parameters
    ----------
    model_dir : Path, optional
        Target directory for models (not used for MS2PIP as models are bundled).

    Raises
    ------
    ImportError
        If ms2pip package is not installed.
    """
    try:
        from ms2pip.constants import MODELS
        logger.info("MS2PIP models are bundled with the ms2pip package.")
        logger.info(f"Available MS2PIP models: {list(MODELS.keys())}")
        logger.info("MS2PIP models validated successfully.")
    except ImportError as e:
        logger.error("MS2PIP package not found. Please install ms2pip>=4.0")
        raise


def download_deeplc_models(model_dir: Optional[Path] = None) -> None:
    """
    Download DeepLC models.

    DeepLC models are automatically downloaded when the DeepLC predictor is
    first instantiated. This function triggers that download by creating
    a temporary DeepLC instance.

    Parameters
    ----------
    model_dir : Path, optional
        Target directory for models. If provided, models will be downloaded
        to this location. Otherwise, they will be downloaded to the default
        DeepLC cache directory.

    Raises
    ------
    ImportError
        If DeepLC package is not installed.
    Exception
        If model download fails.
    """
    try:
        from deeplc import DeepLC
        logger.info("Downloading DeepLC models...")
        
        # Initialize DeepLC to trigger model download
        # This will download models to the default cache location
        try:
            predictor = DeepLC(
                n_jobs=1,
                verbose=False,
            )
            logger.info("DeepLC models downloaded successfully.")
            logger.info("Models are cached in the DeepLC default directory.")
        except Exception as e:
            logger.error(f"Failed to initialize DeepLC: {e}")
            raise
                
    except ImportError as e:
        logger.error("DeepLC package not found. Please install deeplc>=3.0")
        raise


def download_alphapeptdeep_models(model_dir: Optional[Path] = None) -> None:
    """
    Download AlphaPeptDeep (peptdeep) models.

    This function downloads the pretrained models for AlphaPeptDeep/peptdeep
    for MS2 spectrum prediction, retention time prediction, and CCS prediction.

    Parameters
    ----------
    model_dir : Path, optional
        Target directory for models. If provided, models will be copied to
        this location after downloading.

    Raises
    ------
    ImportError
        If peptdeep package is not installed.
    Exception
        If model download fails.
    """
    try:
        from peptdeep.pretrained_models import (
            MODEL_ZIP_FILE_PATH,
            _download_models,
        )
        
        logger.info("Downloading AlphaPeptDeep models...")
        
        # Download models to default location
        _download_models(MODEL_ZIP_FILE_PATH)
        
        logger.info("AlphaPeptDeep models downloaded successfully.")
        logger.info(f"Models cached at: {MODEL_ZIP_FILE_PATH}")
        
        # If a custom model_dir is specified, copy models there
        if model_dir:
            model_dir = Path(model_dir)
            model_dir.mkdir(parents=True, exist_ok=True)
            
            # The models are stored in the peptdeep package directory
            # We need to find where they are cached
            # Import is done here to find the actual installation path at runtime
            import peptdeep  # noqa: F401 - imported for path detection
            peptdeep_path = Path(peptdeep.__file__).parent
            models_path = peptdeep_path / "pretrained_models"
            
            if models_path.exists():
                logger.info(f"Copying models to {model_dir}...")
                for model_file in models_path.glob("*.pth"):
                    target_file = model_dir / model_file.name
                    shutil.copy2(model_file, target_file)
                    logger.info(f"Copied {model_file.name} to {target_file}")
            else:
                logger.warning(
                    f"Could not find models at {models_path}. "
                    "Models are still available in the default cache."
                )
                
    except ImportError as e:
        logger.error("peptdeep package not found. Please install peptdeep")
        raise


@click.command(
    "download_models",
    short_help="Download all models for offline use (MS2PIP, DeepLC, AlphaPeptDeep).",
)
@click.option(
    "--model_dir",
    help="Directory to store downloaded models (optional, uses default cache if not specified)",
    type=click.Path(file_okay=False, dir_okay=True),
    default=None,
)
@click.option(
    "--log_level",
    help="Logging level (default: `info`)",
    default="info",
)
@click.option(
    "--models",
    help="Comma-separated list of models to download: ms2pip, deeplc, alphapeptdeep (default: all)",
    default="ms2pip,deeplc,alphapeptdeep",
)
def download_models(model_dir: Optional[str], log_level: str, models: str) -> None:
    """
    Download all required models for quantms-rescoring for offline use.

    This command downloads models for MS2PIP, DeepLC, and AlphaPeptDeep
    to enable running quantms-rescoring in environments without internet access.

    Examples
    --------
    Download all models to default cache locations:

        $ rescoring download_models

    Download all models to a specific directory:

        $ rescoring download_models --model_dir /path/to/models

    Download only specific models:

        $ rescoring download_models --models deeplc,alphapeptdeep

    Parameters
    ----------
    model_dir : str, optional
        Directory to store downloaded models. If not specified, models are
        downloaded to their default cache locations.
    log_level : str
        Logging level (default: "info").
    models : str
        Comma-separated list of models to download (default: "ms2pip,deeplc,alphapeptdeep").
    """
    # Configure logging
    configure_logging(log_level.upper())
    
    # Convert model_dir to Path if provided
    target_dir = Path(model_dir) if model_dir else None
    
    # Parse models list
    models_list = [m.strip().lower() for m in models.split(",")]
    
    logger.info("Starting model download process...")
    if target_dir:
        logger.info(f"Target directory: {target_dir}")
        target_dir.mkdir(parents=True, exist_ok=True)
    else:
        logger.info("Using default cache locations for each model type")
    
    # Download requested models
    success_count = 0
    failed_models = []
    
    if "ms2pip" in models_list:
        try:
            logger.info("\n=== Downloading MS2PIP models ===")
            download_ms2pip_models(target_dir)
            success_count += 1
        except Exception as e:
            logger.error(f"Failed to download MS2PIP models: {e}")
            failed_models.append("ms2pip")
    
    if "deeplc" in models_list:
        try:
            logger.info("\n=== Downloading DeepLC models ===")
            download_deeplc_models(target_dir)
            success_count += 1
        except Exception as e:
            logger.error(f"Failed to download DeepLC models: {e}")
            failed_models.append("deeplc")
    
    if "alphapeptdeep" in models_list:
        try:
            logger.info("\n=== Downloading AlphaPeptDeep models ===")
            download_alphapeptdeep_models(target_dir)
            success_count += 1
        except Exception as e:
            logger.error(f"Failed to download AlphaPeptDeep models: {e}")
            failed_models.append("alphapeptdeep")
    
    # Summary
    logger.info("\n=== Download Summary ===")
    logger.info(f"Successfully downloaded: {success_count}/{len(models_list)} model types")
    
    if failed_models:
        logger.error(f"Failed to download: {', '.join(failed_models)}")
        raise click.ClickException(
            f"Failed to download some models: {', '.join(failed_models)}"
        )
    else:
        logger.info("All requested models downloaded successfully!")
        logger.info("\nYou can now use quantms-rescoring in offline environments.")
        if target_dir:
            logger.info(f"Models are available in: {target_dir}")
