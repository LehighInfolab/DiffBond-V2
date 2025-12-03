import logging
import os
from datetime import datetime


def setup_logging(verbose=False, log_dir="logs"):
    """
    Configure the root logger to log to both console and a timestamped file.

    After calling this once, you can just use logging.getLogger() anywhere.
    """
    os.makedirs(log_dir, exist_ok=True)

    log_filename = datetime.now().strftime(f"{log_dir}/run_%Y-%m-%d_%H-%M-%S.log")

    # Reset existing handlers to avoid duplicate logs if setup_logging is called again
    logging.getLogger().handlers.clear()

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
    )

    logging.info("Logging initialized. Log file: %s", log_filename)
