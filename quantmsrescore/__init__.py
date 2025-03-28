# Apply warning filter before any imports
from warnings import filterwarnings

filterwarnings(
    "ignore",
    message="OPENMS_DATA_PATH environment variable already exists",
    category=UserWarning,
    module="pyopenms",
)

__version__ = "0.0.6"
