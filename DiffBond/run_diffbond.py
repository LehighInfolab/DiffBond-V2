import os
from pathlib import Path
import logging
import csv
import argparse
from tqdm import tqdm
import subprocess
import concurrent.futures


def setup_logger(log_file="automation.log", level=logging.INFO):
    """
    Create/configure a dedicated logger for the automation script.
    - Writes to ./automation.log
    - Does NOT propagate to root (so it won't get mixed with DiffBond's logging)
    - Safe to call multiple times (won't add duplicate handlers)
    """
    os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)

    logger = logging.getLogger("automation")
    logger.setLevel(level)
    logger.propagate = False  # <- critical: keep it isolated from root/DiffBond

    def _has_logfile_handler(handler):
        base_file = getattr(handler, "baseFilename", None)
        return (
            isinstance(handler, logging.FileHandler)
            and base_file
            and os.path.abspath(base_file) == os.path.abspath(log_file)
        )

    if not any(_has_logfile_handler(h) for h in logger.handlers):
        fh = logging.FileHandler(log_file, mode="a", encoding="utf-8")
        fh.setLevel(level)
        fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s"))
        logger.addHandler(fh)

    return logger


logger = setup_logger()


def make_output_dir(pdb, path):
    try:
        (Path("Results") / path / pdb).mkdir(exist_ok=True)
    except Exception:
        pass


def parse_skempi():
    # Make a skempi directory for current base
    for i in range(8):
        base_dir = Path("results") / f"base-{i}"
        base_dir.mkdir(exist_ok=True)
        logger.info(f"Making new base folder here: {base_dir}")

    # Open skempi csv file and start reading from base1 and split into folders from PDB and subfolders for index
    idx_pdb_dict = {}
    base_idx = 0
    current_pdb = ""

    with open("../../SKEMPI_dataset_download/skempi_v2.csv", "r") as csv_file:
        for index, row in enumerate(csv.reader(csv_file, delimiter=";")):
            if index == 0:  # Skip header
                continue

            pdb = row[0].split("_")[0]
            if pdb != current_pdb:
                current_pdb = pdb
                make_output_dir(pdb, "base-" + str(base_idx))

            if index % 1000 == 0:
                base_idx += 1
                make_output_dir(pdb, "base-" + str(base_idx))

            formatted_index = f"{(index):05d}"
            (Path("results") / f"base-{base_idx}" / pdb / formatted_index).mkdir(exist_ok=True)
            idx_pdb_dict[formatted_index] = f"base-{base_idx}/{pdb}/{formatted_index}"

    return idx_pdb_dict


def calculate_diffbond(file1, file2, output, skempi_csv: str | None = None, row_index: int | None = None):
    # command = ["python", "DiffBond_v2.py", "-i", file1, file2, "-m", "c", "i", "h" ,"s", "-g", "-o", output]
    command = [
        "python",
        "DiffBond_v2.py",
        "-i",
        file1,
        file2,
        "-m",
        "p",
        "-g",
        "--node-level",
        "residue",
        "-o",
        output,
    ]
    if skempi_csv is not None and row_index is not None:
        command.extend(["--skempi-csv", skempi_csv, "--row-index", str(row_index)])
    result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(result.stdout)


def process_folder(folder, path, idx_pdb_dict, is_wt=False):
    if os.path.isfile(path + "/" + folder):
        logger.warning(f"{folder} is a file and not a directory.")
        return

    halfs = [d for d in os.listdir(path + "/" + folder) if d[0] == "H" and len(d) >= 4]

    if len(halfs) != 2:
        logger.error(f"{str(folder)} is missing pdb halfs.")
        return

    halfs.sort(reverse=False)

    try:
        if is_wt:
            pdb = folder.split("/")[-1]
            file1 = f"{path}/{folder}/{halfs[0]}/hydrogensAdded.pdb"
            file2 = f"{path}/{folder}/{halfs[1]}/hydrogensAdded.pdb"
            output = f"wt/{pdb}"
        else:
            file1 = f"{path}/{folder}/{halfs[0]}/half1.pdb"
            file2 = f"{path}/{folder}/{halfs[1]}/half2.pdb"
            output = idx_pdb_dict[folder]

        # Derive row_index from folder name (e.g., '00001' -> 1) for SKEMPI CSV
        row_index = None if is_wt else int(folder)
        skempi_csv = None if is_wt else "datasets/MODIFIED-skempi_v2.csv"

        calculate_diffbond(file1, file2, output, skempi_csv=skempi_csv, row_index=row_index)
    except Exception as error:
        logger.error("Error in calculating diffbond: %s", error)


def _ordered_folders(path, start_from=None):
    """Numeric folders sorted by value; non-numeric left in natural order.
    start_from: int/str like '00325' filters numeric >= that value,
                or a folder *name* to start from exactly.
    """
    try:
        entries = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    except OSError as e:
        logger.error("Could not list directory '%s': %s", path, e, exc_info=True)
        return []

    numeric = sorted([d for d in entries if d.isdigit()], key=lambda s: int(s))
    non_numeric = [d for d in entries if not d.isdigit()]  # keep natural order

    if start_from is None:
        return numeric + non_numeric

    s = str(start_from)
    if s.isdigit():
        start_int = int(s)
        numeric = [d for d in numeric if int(d) >= start_int]
        return numeric + non_numeric
    else:
        ordered = numeric + non_numeric
        if s in ordered:
            return ordered[ordered.index(s) :]
        logger.warning("start_from '%s' not found; processing all.", s)
        return ordered


def diffbond_calc(
    idx_pdb_dict,
    path,
    is_wt=False,
    with_parallel=False,
    start_from=None,
    test_samples_per_bucket=None,
):
    """
    Process subfolders in `path` whose names are zero-padded digits like '00001'.
    Folders are processed in ascending numeric order. Optionally start from a
    particular folder number (int or str), inclusive.

    Args:
        idx_pdb_dict: ...
        path (str): Root directory containing numbered subfolders (e.g., '00001').
        is_wt (bool): Whether to use Results/wt output path.
        with_parallel (bool): (kept for signature compatibility; non-parallel path here)
        start_from (int | str | None): First folder number to process (inclusive).
            Examples: 1, 325, "00001", "00325". If None, starts from the smallest.
        test_samples_per_bucket (int | None): If set, only process this many samples per bucket.
            Useful for testing with a small sample from each range.
    """
    if is_wt:
        try:
            os.makedirs("Results/wt", exist_ok=True)
            logger.info("Making new base folder here: Results/wt")
        except OSError as error:
            logger.error(error)

    folders = _ordered_folders(path, start_from)

    # Limit folders for test mode
    if test_samples_per_bucket is not None and test_samples_per_bucket > 0:
        original_count = len(folders)
        folders = folders[:test_samples_per_bucket]
        logger.info(
            f"TEST MODE: Processing {len(folders)} out of {original_count} folders "
            f"from {'wildtype' if is_wt else path}"
        )
        logger.info(f"Selected folders for testing: {folders}")

    total_folders = len(folders)

    if not with_parallel:
        with tqdm(total=total_folders, desc="Processing folders") as progress_bar:
            for folder in folders:
                try:
                    logger.info(
                        "Processing folder:\n\t%s\n\t%s\n\tis_wt = %s",
                        folder,
                        path,
                        is_wt,
                    )
                    process_folder(folder, path, idx_pdb_dict, is_wt)
                except Exception as exc:
                    logger.error("Error processing %s: %s", folder, exc, exc_info=True)
                finally:
                    progress_bar.update(1)

    else:
        with tqdm(total=total_folders, desc="Processing folders") as progress_bar:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                futures = [executor.submit(process_folder, f, path, idx_pdb_dict, is_wt) for f in folders]
                for future in concurrent.futures.as_completed(futures):
                    try:
                        future.result()
                    except Exception as exc:
                        logger.error("Generated an exception: %s", exc)
                    finally:
                        progress_bar.update(1)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run DiffBond calculations on SKEMPI dataset with optional test mode.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--test",
        type=int,
        metavar="N",
        help="Test mode: process only N samples from each bucket (wildtype and each base range). "
        "Useful for quick testing. Default: process all samples.",
    )

    parser.add_argument(
        "--base-range",
        type=int,
        nargs=2,
        metavar=("START", "END"),
        help="Process only base folders in the range [START, END) (0-indexed). "
        "Example: --base-range 0 2 processes base-0 and base-1. Default: process all bases.",
    )

    parser.add_argument(
        "--skip-wt",
        action="store_true",
        help="Skip wildtype processing.",
    )

    parser.add_argument(
        "--wt-path",
        type=str,
        default="../../SKEMPI_dataset_download/wt",
        help="Path to wildtype directory. Default: ../../SKEMPI_dataset_download/wt",
    )

    parser.add_argument(
        "--base-path-template",
        type=str,
        default="../../MT_Processing_Archive/base-{i}",
        help="Template for base folder paths. Use {i} as placeholder for base number. "
        "Default: ../../MT_Processing_Archive/base-{i}",
    )

    return parser.parse_args()


def main():
    args = parse_arguments()
    reader = parse_skempi()

    test_samples = args.test
    if test_samples:
        logger.info("=" * 60 + "\n" f"TEST MODE ENABLED: Processing {test_samples} samples per bucket\n" "=" * 60)

    # Process wild type (if not skipped)
    if not args.skip_wt:
        wt_path = args.wt_path
        logger.info(
            "---------------------------------------------------------------\n"
            f"Starting diffbond_calc on wildtype on path: {wt_path}\n"
            "--------------------------------------------------------------- "
        )
        diffbond_calc(reader, wt_path, is_wt=True, test_samples_per_bucket=test_samples)
        logger.info(
            "---------------------------------------------------------------\n"
            f"Finished diffbond_calc on wildtype on path: {wt_path}\n"
            "--------------------------------------------------------------- "
        )
    else:
        logger.info("Skipping wildtype processing (--skip-wt flag set)")

    # Determine base range to process
    if args.base_range:
        base_start, base_end = args.base_range
        if base_start < 0 or base_end <= base_start:
            logger.error(f"Invalid base range: [{base_start}, {base_end})")
            return
        base_indices = range(base_start, base_end)
        logger.info(f"Processing bases in range: {base_start} to {base_end - 1}")
    else:
        # Default: process all bases (0-7)
        base_indices = range(8)

    # Process base folders
    for i in base_indices:
        base_path = args.base_path_template.format(i=i)
        logger.info(
            "---------------------------------------------------------------\n"
            f"Starting diffbond_calc on base_path: {base_path}\n"
            "--------------------------------------------------------------- "
        )
        diffbond_calc(reader, base_path, test_samples_per_bucket=test_samples)
        logger.info(
            "---------------------------------------------------------------\n"
            f"Finished diffbond_calc on base_path: {base_path}\n"
            "--------------------------------------------------------------- "
        )


if __name__ == "__main__":
    main()

# =====================================================================
# EXAMPLE COMMANDS
# =====================================================================

# Process all data (default: wildtype + all 8 bases)
# python run_diffbond.py

# Process all bases but skip wildtype
# python run_diffbond.py --skip-wt

# Test mode: process 3 samples from each bucket (wildtype + base-0 through base-7)
# python run_diffbond.py --test 3

# Test mode skipping wildtype rtimojuhb b jgnvm n bm    im gay x100

# python run_diffbond.py --test 10 --skip-wt

# Process only the first 1000 samples (base-0 only, rows 1-1000)
# python run_diffbond.py --base-range 0 1 --skip-wt

# Process first 2000 samples (base-0 and base-1, rows 1-2000)
# python run_diffbond.py --base-range 0 2 --skip-wt

# Process a specific range: base-2 through base-5 (rows 2001-5000)
# python run_diffbond.py --base-range 2 6 --skip-wt

# Process only wildtype (skip all mutation bases)
# Note: Currently there's no way to process ONLY wildtype without also processing bases.
# The --skip-wt flag skips wildtype, not the bases.
