import logging
import os
import sys
import io
from datetime import datetime


def setup_logging(verbose=False, create_file_handlers=True):
    """
    Configure the root logger to log to both console and a file.
    Also sets up a separate failure logger for edge calculation failures.

    Args:
        verbose: Enable verbose (DEBUG) logging
        create_file_handlers: If False, only set up console logging and UTF-8 handling,
            without creating file handlers. Useful for initial setup before log file paths are known.

    After calling this once, you can just use logging.getLogger() anywhere.
    Use logging.getLogger("failures") for failure-specific logging.

    Log files are always written to:
    - Main log: logs/diffbond_all.log
    - Failure log: logs/diffbond_failures.log

    Note: Multiple processes can safely call this - they will all append to the same files.
    """
    root_logger = logging.getLogger()
    failure_logger = logging.getLogger("failures")

    # Fixed log file paths
    log_dir = "logs"
    log_filename = os.path.join(log_dir, "diffbond_all.log")
    failure_log_filename = os.path.join(log_dir, "diffbond_failures.log")

    # If create_file_handlers is False, skip file handler creation (for initial setup)
    if not create_file_handlers:
        # Initial setup - only configure console logging and UTF-8, no file handlers
        # Reconfigure std streams to use UTF-8 with replacement for errors
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
            sys.stderr.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

        # Set up console handler if not already present
        has_stream_handler = any(
            isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
            for h in root_logger.handlers
        )
        if not has_stream_handler:
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
            root_logger.addHandler(stream_handler)
            root_logger.setLevel(logging.DEBUG if verbose else logging.INFO)
        return

    # Normal setup - create file handlers
    # Ensure directory exists
    os.makedirs(log_dir, exist_ok=True)

    # Get loggers (already retrieved above, but keep for clarity)

    def _has_file_handler(logger_instance, target_file):
        """Check if logger already has a FileHandler for the target file."""
        abs_target = os.path.abspath(target_file)
        for handler in logger_instance.handlers:
            if isinstance(handler, logging.FileHandler):
                abs_base = os.path.abspath(handler.baseFilename)
                if abs_base == abs_target:
                    return True
        return False

    # Check if we need to reconfigure handlers
    needs_root_setup = not _has_file_handler(root_logger, log_filename)
    needs_failure_setup = not _has_file_handler(failure_logger, failure_log_filename)

    if needs_root_setup:
        # Remove existing file handlers if we're reconfiguring to a different file
        # Keep stream handlers intact
        for handler in list(root_logger.handlers):
            if isinstance(handler, logging.FileHandler):
                root_logger.removeHandler(handler)
                handler.close()

        # Ensure file handler writes UTF-8 so logs are portable and can contain unicode
        # Use 'a' mode to append if file exists
        file_handler = logging.FileHandler(log_filename, encoding="utf-8", mode="a")
        file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))

        # Reconfigure std streams to use UTF-8 with replacement for errors so
        # non-encodable characters don't raise UnicodeEncodeError on Windows consoles
        try:
            # Python 3.7+: TextIOBase.reconfigure is available
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
            sys.stderr.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            # Fallback: wrap the buffer in a TextIOWrapper
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

        # Check if we already have a stream handler
        has_stream_handler = any(
            isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
            for h in root_logger.handlers
        )
        if not has_stream_handler:
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
            root_logger.addHandler(stream_handler)

        # Add the file handler
        root_logger.addHandler(file_handler)
        root_logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    else:
        # Just update the log level if handlers already exist for the correct file
        root_logger.setLevel(logging.DEBUG if verbose else logging.INFO)

    if needs_failure_setup:
        # Set up separate failure logger
        failure_logger.setLevel(logging.WARNING)  # Only warnings and errors
        failure_logger.propagate = False  # Don't propagate to root logger

        failure_file_handler = logging.FileHandler(failure_log_filename, encoding="utf-8", mode="a")
        failure_file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        failure_logger.addHandler(failure_file_handler)

    if needs_root_setup or needs_failure_setup:
        logging.info("Logging initialized. Main log file: %s", log_filename)
        logging.info("Failure log file: %s", failure_log_filename)
