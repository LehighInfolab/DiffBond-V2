"""
Run DiffBond calculations on SKEMPI dataset.

Command-line Parameters:
    --test N
        Process only N samples from each bucket (wildtype and each base range).
        Example: --test 3

    --base-range START END
        Process only base folders in the range [START, END) (0-indexed).
        Example: --base-range 0 2 processes base-0 and base-1.
        Cannot be used together with --row-range.

    --row-range START END
        Process only rows in the range [START, END) (exclusive end).
        Example: --row-range 4739 5000 processes rows 4739-4999.
        Cannot be used together with --base-range.

    --skip-wt
        Skip wildtype processing entirely.

    --skip-mutant
        Skip mutant distance calculation even when SKEMPI CSV is available.

    --wt-pdb PDB [PDB ...]
        Process only these wildtype PDB code(s).
        Example: --wt-pdb 3UIH 4G0N or --wt-pdb 3UIH,4G0N

Example Commands:
    # Process all data (default: wildtype + all 8 bases)
    python run_diffbond.py

    # Process all bases but skip wildtype
    python run_diffbond.py --skip-wt

    # Test mode: process 3 samples from each bucket
    python run_diffbond.py --test 3

    # Process only the first 1000 samples (base-0 only)
    python run_diffbond.py --base-range 0 1 --skip-wt

    # Process exact row range: rows 4739 to 4999
    python run_diffbond.py --row-range 4739 5000 --skip-wt

    # Process only specific wildtype PDBs
    python run_diffbond.py --wt-pdb 3UIH 4G0N

    # Skip mutant distance calculation
    python run_diffbond.py --skip-mutant
"""

import os
from pathlib import Path
import logging
import csv
import argparse
from tqdm import tqdm
import subprocess


def setup_logger(log_file="logs/run_diffbond.log", level=logging.INFO):
    """Create/configure a dedicated logger for the automation script."""
    os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)

    logger = logging.getLogger("automation")
    logger.setLevel(level)
    logger.propagate = False

    def _has_logfile_handler(handler):
        base_file = getattr(handler, "baseFilename", None)
        return isinstance(handler, logging.FileHandler) and (
            base_file and os.path.abspath(base_file) == os.path.abspath(log_file)
        )

    if not any(_has_logfile_handler(h) for h in logger.handlers):
        fh = logging.FileHandler(log_file, mode="a", encoding="utf-8")
        fh.setLevel(level)
        fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s"))
        logger.addHandler(fh)

    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        ch = logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(ch)

    return logger


logger = setup_logger()


def parse_skempi():
    """Parse SKEMPI CSV and create output directory structure."""
    for i in range(8):
        (Path("Results") / f"base-{i}").mkdir(parents=True, exist_ok=True)

    idx_pdb_dict = {}
    base_idx = 0
    current_pdb = ""

    with open("../../SKEMPI_dataset_download/skempi_v2.csv", "r") as csv_file:
        for index, row in enumerate(csv.reader(csv_file, delimiter=";")):
            if index == 0:
                continue

            pdb = row[0].split("_")[0]
            if pdb != current_pdb:
                current_pdb = pdb
                (Path("Results") / f"base-{base_idx}" / pdb).mkdir(parents=True, exist_ok=True)

            if index % 1000 == 0:
                base_idx += 1
                (Path("Results") / f"base-{base_idx}" / pdb).mkdir(parents=True, exist_ok=True)

            formatted_index = f"{index:05d}"
            (Path("Results") / f"base-{base_idx}" / pdb / formatted_index).mkdir(parents=True, exist_ok=True)
            idx_pdb_dict[formatted_index] = f"base-{base_idx}/{pdb}/{formatted_index}"

    return idx_pdb_dict


def calculate_diffbond(file1, file2, output, skempi_csv=None, row_index=None, skip_mutant=False):
    """Run DiffBond_v2.py as a subprocess."""
    command = [
        "python",
        "DiffBond_v2.py",
        "-i",
        file1,
        file2,
        "-m",
        "c,i,h,s,p,b,m,n",
        "-g",
        "--node-level",
        "residue",
        "-o",
        output,
    ]
    if not skip_mutant and skempi_csv and row_index is not None:
        command.extend(["--skempi-csv", skempi_csv, "--row-index", str(row_index)])

    try:
        result = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
        )
        if result.stdout:
            print(result.stdout)
    except subprocess.CalledProcessError as e:
        error_msg = f"DiffBond_v2.py failed for {output}"
        if e.stdout:
            error_msg += f"\nSubprocess output:\n{e.stdout}"
        logger.error(error_msg)
        raise


def process_folder(folder, path, idx_pdb_dict, is_wt=False, skip_mutant=False):
    """Process a single folder by running DiffBond calculations."""
    folder_path = Path(path) / folder
    if not folder_path.is_dir():
        logger.warning(f"{folder} is not a directory.")
        return

    halfs = sorted([d for d in os.listdir(folder_path) if d.startswith("H") and len(d) >= 4])
    if len(halfs) != 2:
        logger.error(f"{folder} is missing PDB halves.")
        return

    if is_wt:
        pdb = folder.split("/")[-1]
        file1 = folder_path / halfs[0] / "hydrogensAdded.pdb"
        file2 = folder_path / halfs[1] / "hydrogensAdded.pdb"
        output = f"wt/{pdb}"
        skempi_csv = None
        row_index = None
    else:
        file1 = folder_path / halfs[0] / "half1.pdb"
        file2 = folder_path / halfs[1] / "half2.pdb"
        output = idx_pdb_dict[folder]
        skempi_csv = "datasets/MODIFIED-skempi_v2.csv"
        row_index = int(folder)

    calculate_diffbond(
        str(file1),
        str(file2),
        output,
        skempi_csv=skempi_csv,
        row_index=row_index,
        skip_mutant=skip_mutant,
    )


def _ordered_folders(path, start_from=None, end_at=None):
    """Get folders sorted numerically, optionally filtered by range."""
    try:
        entries = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    except OSError as e:
        logger.error("Could not list directory '%s': %s", path, e, exc_info=True)
        return []

    numeric = sorted([d for d in entries if d.isdigit()], key=int)
    non_numeric = sorted([d for d in entries if not d.isdigit()])
    ordered = numeric + non_numeric

    if start_from is not None:
        start_int = int(start_from) if str(start_from).isdigit() else None
        if start_int is not None:
            ordered = [d for d in ordered if not d.isdigit() or int(d) >= start_int]
        elif str(start_from) in ordered:
            ordered = ordered[ordered.index(str(start_from)) :]

    if end_at is not None and str(end_at).isdigit():
        end_int = int(end_at)
        ordered = [d for d in ordered if not d.isdigit() or int(d) < end_int]

    return ordered


def diffbond_calc(
    idx_pdb_dict,
    path,
    is_wt=False,
    start_from=None,
    end_at=None,
    test_samples_per_bucket=None,
    skip_mutant=False,
    wt_pdbs=None,
):
    """Process subfolders in path, running DiffBond calculations."""
    if is_wt:
        os.makedirs("Results/wt", exist_ok=True)

    folders = _ordered_folders(path, start_from, end_at)

    if is_wt and wt_pdbs:
        wanted = {p.upper() for p in wt_pdbs if p}
        folders = [f for f in folders if f.upper() in wanted]
        missing = sorted(wanted - {f.upper() for f in folders})
        if missing:
            logger.warning("Requested WT PDB(s) not found: %s", ", ".join(missing))
        logger.info("Wildtype filter: %d folder(s) selected", len(folders))

    if test_samples_per_bucket:
        original_count = len(folders)
        folders = folders[:test_samples_per_bucket]
        logger.info(f"TEST MODE: Processing {len(folders)} out of {original_count} folders")

    with tqdm(total=len(folders), desc="Processing folders") as progress_bar:
        for folder in folders:
            try:
                process_folder(folder, path, idx_pdb_dict, is_wt, skip_mutant=skip_mutant)
            except Exception as exc:
                logger.error("Error processing %s: %s", folder, exc, exc_info=True)
            finally:
                progress_bar.update(1)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run DiffBond calculations on SKEMPI dataset.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--test",
        type=int,
        metavar="N",
        help="Test mode: process only N samples from each bucket.",
    )

    parser.add_argument(
        "--base-range",
        type=int,
        nargs=2,
        metavar=("START", "END"),
        help="Process only base folders in the range [START, END) (0-indexed).",
    )

    parser.add_argument(
        "--row-range",
        type=int,
        nargs=2,
        metavar=("START", "END"),
        help="Process only rows in the range [START, END) (exclusive end).",
    )

    parser.add_argument("--skip-wt", action="store_true", help="Skip wildtype processing.")

    parser.add_argument(
        "--wt-pdb",
        type=str,
        nargs="+",
        metavar="PDB",
        help="Process only these wildtype PDB code(s).",
    )

    parser.add_argument(
        "--skip-mutant",
        action="store_true",
        help="Skip mutant distance calculation.",
    )

    return parser.parse_args()


def main():
    args = parse_arguments()

    WT_PATH = "../../SKEMPI_dataset_download/WT_PPI_processing"
    BASE_PATH_TEMPLATE = "../../MT_Processing_Archive/base-{i}"

    reader = parse_skempi()

    wt_pdbs = None
    if args.wt_pdb:
        wt_pdbs = []
        for token in args.wt_pdb:
            wt_pdbs.extend([p.strip() for p in str(token).split(",") if p.strip()])

    os.makedirs("logs", exist_ok=True)
    logger.info("DiffBond instances will log to: logs/diffbond_all.log")

    if args.test:
        logger.info(f"TEST MODE ENABLED: Processing {args.test} samples per bucket")

    if args.base_range and args.row_range:
        logger.error("Cannot use --base-range and --row-range together.")
        return

    row_start, row_end = args.row_range if args.row_range else (None, None)
    if args.row_range:
        if row_start < 0 or row_end <= row_start:
            logger.error(f"Invalid row range: [{row_start}, {row_end})")
            return
        logger.info(f"Processing rows in range: {row_start} to {row_end - 1}")

    if not args.skip_wt:
        logger.info(f"Starting diffbond_calc on wildtype: {WT_PATH}")
        diffbond_calc(
            reader,
            WT_PATH,
            is_wt=True,
            start_from=row_start,
            end_at=row_end,
            test_samples_per_bucket=args.test,
            skip_mutant=args.skip_mutant,
            wt_pdbs=wt_pdbs,
        )
        logger.info("Finished diffbond_calc on wildtype")

    if args.base_range:
        if args.base_range[0] < 0 or args.base_range[1] <= args.base_range[0]:
            logger.error(f"Invalid base range: [{args.base_range[0]}, {args.base_range[1]})")
            return
        base_indices = range(*args.base_range)
        logger.info(f"Processing bases in range: {args.base_range[0]} to {args.base_range[1] - 1}")
    else:
        base_indices = range(8)

    for i in base_indices:
        base_path = BASE_PATH_TEMPLATE.format(i=i)
        logger.info(f"Starting diffbond_calc on base_path: {base_path}")
        diffbond_calc(
            reader,
            base_path,
            start_from=row_start,
            end_at=row_end,
            test_samples_per_bucket=args.test,
            skip_mutant=args.skip_mutant,
        )
        logger.info(f"Finished diffbond_calc on base_path: {base_path}")


if __name__ == "__main__":
    main()
