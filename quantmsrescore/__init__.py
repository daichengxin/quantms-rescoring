# Handle environment variables and warnings before any imports
import os
import sys
from warnings import filterwarnings

# Store the original stdout to restore it later
original_stdout = sys.stdout
# Temporarily redirect stdout to suppress any output during import
sys.stdout = open(os.devnull, 'w')

# Apply multiple warning filters with different patterns to catch all variations
filterwarnings("ignore", message=".*OPENMS_DATA_PATH.*", category=UserWarning)
filterwarnings("ignore", message=".*OPENMS_DATA_PATH.*")
filterwarnings("ignore", category=UserWarning, module="pyopenms")
filterwarnings("ignore", module="pyopenms")

# Try to handle the environment variable directly
# If we're getting a warning about it already existing, we can temporarily unset it
# and let pyOpenMS set it to its default, then restore it after import
openms_data_path = os.environ.get("OPENMS_DATA_PATH")
if openms_data_path:
    # Temporarily unset the variable
    del os.environ["OPENMS_DATA_PATH"]

# After all our preparations, we can now safely import pyOpenMS
try:
    import pyopenms
except ImportError:
    # If import fails, it's not critical at this point
    pass

# Restore the environment variable if we unset it
if openms_data_path:
    os.environ["OPENMS_DATA_PATH"] = openms_data_path

# Restore stdout
sys.stdout.close()
sys.stdout = original_stdout

__version__ = "0.0.6"
