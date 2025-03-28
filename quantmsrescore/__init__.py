import warnings
import sys
import atexit

# Capture and filter warnings at the lowest level possible
original_showwarning = warnings.showwarning

def custom_showwarning(message, category, filename, lineno, file=None, line=None):
    msg_str = str(message)
    # Filter out specific warnings
    if any(pattern in msg_str for pattern in [
        "OPENMS_DATA_PATH",
        "Could not add the following value",
        "Skipping the following (not in library)",
        "DeepLC tried to set intra op threads"
    ]):
        return  # Suppress the warning
    # For all other warnings, use the original handler
    return original_showwarning(message, category, filename, lineno, file, line)

warnings.showwarning = custom_showwarning

# Also set up standard warning filters
warnings.filterwarnings("ignore", message=".*OPENMS_DATA_PATH.*")
warnings.filterwarnings("ignore", message=".*Could not add the following value.*")
warnings.filterwarnings("ignore", message=".*Skipping the following \(not in library\).*")
warnings.filterwarnings("ignore", message=".*DeepLC tried to set intra op threads.*")

# Create a class to filter stderr output
class StderrFilter:
    def __init__(self):
        self.original_stderr = sys.stderr

    def __enter__(self):
        sys.stderr = self
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stderr = self.original_stderr

    def write(self, data):
        # Filter out specific warning messages
        if any(pattern in data for pattern in [
            "OPENMS_DATA_PATH",
            "Could not add the following value",
            "Skipping the following (not in library)",
            "DeepLC tried to set intra op threads"
        ]):
            return
        # Write other messages to the original stderr
        self.original_stderr.write(data)

    def flush(self):
        self.original_stderr.flush()

    def isatty(self):
        return hasattr(self.original_stderr, 'isatty') and self.original_stderr.isatty()

    def fileno(self):
        return self.original_stderr.fileno()

# Set up the stderr filter
_stderr_filter = StderrFilter()
_stderr_filter.__enter__()

# Register cleanup function to restore stderr
def _cleanup_stderr():
    _stderr_filter.__exit__(None, None, None)

atexit.register(_cleanup_stderr)

__version__ = "0.0.6"