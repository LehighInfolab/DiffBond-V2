"""
Rerun DiffBond calculations for mutations that failed due to parsing issues.

This script:
1. Reads failure log to extract failed mutation tokens
2. Finds corresponding row indices in SKEMPI CSV
3. Reruns DiffBond calculations for those specific folders
4. Uses the same logging system as run_diffbond.py

Command-line Parameters:
    --failure-log PATH
        Path to failure log file (default: logs/diffbond_failures.log)

    --csv PATH
        Path to SKEMPI CSV file (default: datasets/MODIFIED-skempi_v2.csv)

Example Commands:
    # Rerun failed mutations using default paths
    python rerun_failed_mutations.py

    # Use custom failure log
    python rerun_failed_mutations.py --failure-log logs/custom_failures.log
"""

import os
import re
import csv
import argparse
from pathlib import Path
from tqdm import tqdm

from run_diffbond import calculate_diffbond, parse_skempi, setup_logger

logger = setup_logger(log_file="logs/rerun_failed_mutations.log")


def extract_failed_mutations(failure_log_path: str) -> set:
    """Extract unique failed mutation tokens from failure log."""
    mutations = set()
    with open(failure_log_path, "r") as f:
        for line in f:
            match = re.search(r"Failed to parse mutation token '([^']+)'", line)
            if match:
                mutations.add(match.group(1))
    return mutations


def find_row_indices_for_mutations(csv_path: str, mutations: set) -> list:
    """Find row indices in SKEMPI CSV that contain the given mutations."""
    indices = []
    with open(csv_path, "r") as f:
        reader = csv.reader(f, delimiter=";")
        for row_idx, row in enumerate(reader):
            if len(row) > 1 and row[1].strip():
                row_mutations = [m.strip() for m in row[1].strip().split(",")]
                for mut in mutations:
                    if mut in row_mutations:
                        folder_index = row_idx + 1
                        indices.append(folder_index)
                        break
    return sorted(set(indices))


def rerun_failed_mutations(
    failure_log_path: str = "logs/diffbond_failures.log",
    csv_path: str = "datasets/MODIFIED-skempi_v2.csv",
) -> None:
    """Rerun DiffBond calculations for mutations that failed due to parsing issues."""
    BASE_PATH_TEMPLATE = "../../MT_Processing_Archive/base-{i}"

    logger.info("=" * 80)
    logger.info("RERUN FAILED MUTATIONS: Starting...")
    logger.info("=" * 80)

    if not os.path.exists(failure_log_path):
        logger.error(f"Failure log not found: {failure_log_path}")
        return

    logger.info(f"Reading failure log from: {failure_log_path}")
    mutations = extract_failed_mutations(failure_log_path)
    logger.info(f"Found {len(mutations)} unique failed mutations")

    if not mutations:
        logger.info("No failed mutations found. Nothing to rerun.")
        return

    if not os.path.exists(csv_path):
        logger.error(f"CSV file not found: {csv_path}")
        return

    logger.info(f"Searching for mutations in CSV: {csv_path}")
    unique_indices = find_row_indices_for_mutations(csv_path, mutations)
    logger.info(f"Found {len(unique_indices)} unique row indices to reprocess")

    if len(unique_indices) <= 20:
        logger.info(f"Indices: {unique_indices}")
    else:
        logger.info(f"First 20 indices: {unique_indices[:20]}... (total: {len(unique_indices)})")

    logger.info("Parsing SKEMPI structure...")
    idx_pdb_dict = parse_skempi()
    logger.info(f"Parsed {len(idx_pdb_dict)} folder mappings")

    os.makedirs("logs", exist_ok=True)
    logger.info("DiffBond instances will log to: logs/diffbond_all.log")
    logger.info("DiffBond failures will log to: logs/diffbond_failures.log")

    indices_by_base = {}
    for idx in unique_indices:
        base_num = (idx - 1) // 1000
        if base_num not in indices_by_base:
            indices_by_base[base_num] = []
        indices_by_base[base_num].append(idx)

    logger.info(f"Processing indices across {len(indices_by_base)} base folders")
    for base_num in sorted(indices_by_base.keys()):
        logger.info(f"  base-{base_num}: {len(indices_by_base[base_num])} indices")

    for base_num, base_indices in indices_by_base.items():
        base_path = Path(BASE_PATH_TEMPLATE.format(i=base_num))
        if not base_path.exists():
            logger.warning(f"Base path does not exist: {base_path}, skipping")
            continue

        logger.info("=" * 80)
        logger.info(f"Processing base-{base_num} with {len(base_indices)} indices")
        logger.info(f"Base path: {base_path}")
        logger.info("=" * 80)

        all_folders = [d for d in os.listdir(base_path) if os.path.isdir(base_path / d) and d.isdigit()]
        target_folders = sorted([f for f in all_folders if int(f) in base_indices], key=int)

        logger.info(f"Found {len(target_folders)} matching folders in base-{base_num}")

        if not target_folders:
            logger.warning(f"No matching folders found in base-{base_num}, skipping")
            continue

        for folder in tqdm(target_folders, desc=f"Processing base-{base_num}"):
            try:
                folder_path = base_path / folder
                halfs = sorted([d for d in os.listdir(folder_path) if d.startswith("H") and len(d) >= 4])
                if len(halfs) != 2:
                    logger.error(f"{folder} is missing PDB halves.")
                    continue

                file1 = folder_path / halfs[0] / "half1.pdb"
                file2 = folder_path / halfs[1] / "half2.pdb"
                output = idx_pdb_dict[folder]
                row_index = int(folder)

                logger.info(f"Rerunning: base-{base_num}/{folder} -> {output}")
                calculate_diffbond(
                    str(file1),
                    str(file2),
                    output,
                    skempi_csv=csv_path,
                    row_index=row_index,
                    skip_mutant=False,
                )
                logger.info(f"Successfully processed folder {folder} in base-{base_num}")
            except Exception as e:
                logger.error(f"Error processing folder {folder} in base-{base_num}: {e}", exc_info=True)

    logger.info("=" * 80)
    logger.info("Finished reprocessing failed mutations")
    logger.info("=" * 80)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Rerun DiffBond calculations for mutations that failed due to parsing issues.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--failure-log",
        type=str,
        default="logs/diffbond_failures.log",
        help="Path to failure log file",
    )

    parser.add_argument(
        "--csv",
        type=str,
        default="datasets/MODIFIED-skempi_v2.csv",
        help="Path to SKEMPI CSV file",
    )

    return parser.parse_args()


def main():
    """Main function for rerun-failed-mutations mode."""
    args = parse_arguments()

    os.makedirs("logs", exist_ok=True)
    logger.info("DiffBond instances will log to: logs/diffbond_all.log")
    logger.info("DiffBond failures will log to: logs/diffbond_failures.log")

    rerun_failed_mutations(failure_log_path=args.failure_log, csv_path=args.csv)


if __name__ == "__main__":
    main()
