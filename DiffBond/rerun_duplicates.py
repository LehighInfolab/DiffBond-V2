"""
Detect duplicate nodes in edge files and rerun DiffBond calculations.

This script:
1. Detects folders with duplicate nodes (same chain and resSeq but different resName)
2. Reruns DiffBond calculations for those folders
3. Uses the same logging system as run_diffbond.py

Command-line Parameters:
    --results-root PATH
        Root Results directory for duplicate detection (default: ./Results)

    --skip-mutant
        Skip mutant distance calculation even when SKEMPI CSV is available

    --print-only
        Only print found duplicates without rerunning calculations

Example Commands:
    # Detect and rerun folders with duplicates
    python rerun_duplicates.py

    # Only print duplicates without rerunning
    python rerun_duplicates.py --print-only

    # Skip mutant distance calculation
    python rerun_duplicates.py --skip-mutant
"""

import os
import csv
import ast
import argparse
from pathlib import Path
from typing import Optional, List, Tuple, Dict
from collections import defaultdict
from tqdm import tqdm

from run_diffbond import calculate_diffbond, parse_skempi, setup_logger

logger = setup_logger()


def parse_node_string(node_str: str) -> Optional[Tuple[str, str, str]]:
    """Parse a node string representation to extract (chain, resSeq, resName)."""
    if not node_str:
        return None

    node_str = node_str.strip()
    if node_str.startswith('"') and node_str.endswith('"'):
        node_str = node_str[1:-1]
    elif node_str.startswith("'") and node_str.endswith("'"):
        node_str = node_str[1:-1]

    try:
        node_tuple = ast.literal_eval(node_str)
        if isinstance(node_tuple, tuple) and len(node_tuple) >= 3:
            return (str(node_tuple[0]).strip(), str(node_tuple[1]).strip(), str(node_tuple[2]).strip())
    except (ValueError, SyntaxError):
        pass
    return None


def read_edge_file(edge_file: Path) -> Tuple[Dict[int, Tuple[str, str, str]], List[List[int]], bool]:
    """Read an edge file and return nodes and edges."""
    node_index = {}
    edges = []
    has_distance = False

    if not edge_file.exists():
        return node_index, edges, has_distance

    with open(edge_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        rows = list(reader)

    if not rows:
        return node_index, edges, has_distance

    edges_start_idx = None
    for i, row in enumerate(rows):
        if len(row) >= 2 and row[0].strip().lower() == "source" and row[1].strip().lower() == "target":
            edges_start_idx = i
            has_distance = len(row) >= 3 and row[2].strip().lower() == "distance"
            break

    if edges_start_idx is None:
        edges_start_idx = len(rows)

    for i in range(1, edges_start_idx):
        if len(rows[i]) < 2:
            continue
        try:
            index = int(rows[i][0].strip())
            node_tuple = parse_node_string(rows[i][1].strip())
            if node_tuple:
                node_index[index] = node_tuple
        except (ValueError, IndexError):
            continue

    for i in range(edges_start_idx + 1, len(rows)):
        if len(rows[i]) < 2:
            continue
        try:
            source = int(rows[i][0].strip())
            target = int(rows[i][1].strip())
            if has_distance and len(rows[i]) >= 3:
                edges.append([source, target, float(rows[i][2].strip())])
            else:
                edges.append([source, target])
        except (ValueError, IndexError):
            continue

    return node_index, edges, has_distance


def find_duplicate_nodes(node_index: Dict[int, Tuple[str, str, str]]) -> Dict[Tuple[str, str], List[int]]:
    """Find duplicate nodes where same (chain, resSeq) appears with different resName."""
    position_to_indices = defaultdict(list)

    for index, (chain, resseq, resname) in node_index.items():
        position_key = (chain, resseq)
        position_to_indices[position_key].append((index, resname))

    duplicates = {}
    for pos, index_resname_list in position_to_indices.items():
        if len(index_resname_list) > 1:
            resnames = set(resname for _, resname in index_resname_list)
            if len(resnames) > 1:
                duplicates[pos] = [index for index, _ in index_resname_list]

    return duplicates


def detect_duplicates_in_file(edge_file: Path) -> bool:
    """Detect if an edge file has duplicate nodes."""
    node_index, _, _ = read_edge_file(edge_file)
    if not node_index:
        return False

    duplicates = find_duplicate_nodes(node_index)
    return len(duplicates) > 0


def get_duplicates_from_file(edge_file: Path) -> Dict[Tuple[str, str], List[int]]:
    """Get duplicate nodes from an edge file."""
    node_index, _, _ = read_edge_file(edge_file)
    if not node_index:
        return {}

    return find_duplicate_nodes(node_index)


def get_row_index_from_path(results_path: Path, edge_file_dir: Path) -> Optional[int]:
    """Extract row index from directory path."""
    rel_path = edge_file_dir.relative_to(results_path)
    parts = rel_path.parts
    if len(parts) >= 3:
        return int(parts[-1])
    return None


def find_folders_with_duplicates(
    results_root: Optional[Path] = None,
    test_folder: Optional[Path] = None,
) -> List[Tuple[Path, Optional[int], Optional[Path]]]:
    """Find all folders containing edge files with duplicate nodes."""
    folders_with_duplicates = []

    if test_folder is not None:
        results_path = Path(test_folder)
        if not results_path.exists():
            logger.error(f"Test folder not found: {results_path}")
            return folders_with_duplicates
        logger.info(f"[TEST MODE] Checking edge files in: {results_path.absolute()}")
    else:
        if results_root is None:
            results_path = Path("Results")
        else:
            results_path = Path(results_root)

        if not results_path.exists():
            logger.error(f"Results directory not found: {results_path}")
            return folders_with_duplicates

        logger.info(f"Checking edge files in: {results_path.absolute()}")

    if test_folder:
        edge_files = list(results_path.glob("*_edges.txt"))
    else:
        edge_files = list(results_path.rglob("*_edges.txt"))

    filtered_edge_files = []
    for edge_file in edge_files:
        if edge_file.name == "mutant_edges.txt":
            rel_path = edge_file.relative_to(results_path)
            parts = rel_path.parts
            if len(parts) >= 3 and parts[0] == "wt" and parts[-1] == "mutant_edges.txt":
                continue
        filtered_edge_files.append(edge_file)

    if not filtered_edge_files:
        logger.info("No *_edges.txt files found.")
        return folders_with_duplicates

    logger.info(f"Found {len(filtered_edge_files)} edge file(s) to check...")

    files_by_dir = defaultdict(list)
    for edge_file in filtered_edge_files:
        files_by_dir[edge_file.parent].append(edge_file)

    dirs_with_duplicates = set()

    for edge_file_dir, files in tqdm(files_by_dir.items(), desc="Checking for duplicates"):
        contact_file = edge_file_dir / "contact_edges.txt"
        if contact_file in files and detect_duplicates_in_file(contact_file):
            dirs_with_duplicates.add(edge_file_dir)
            continue

        for edge_file in files:
            if edge_file.name != "contact_edges.txt" and detect_duplicates_in_file(edge_file):
                dirs_with_duplicates.add(edge_file_dir)
                break

    for edge_file_dir in dirs_with_duplicates:
        row_index = (
            get_row_index_from_path(results_path, edge_file_dir) if results_path in edge_file_dir.parents else None
        )
        folders_with_duplicates.append((edge_file_dir, row_index, results_path))

    return folders_with_duplicates


def rerun_folders_with_duplicates(
    folders_with_duplicates: List[Tuple[Path, Optional[int], Path]],
    skip_mutant: bool = False,
) -> None:
    """Rerun DiffBond calculations for folders that have duplicate nodes in edge files."""
    if not folders_with_duplicates:
        logger.info("No folders with duplicates found. Nothing to rerun.")
        return

    BASE_PATH_TEMPLATE = "../../MT_Processing_Archive/base-{i}"
    WT_PATH = "../../SKEMPI_dataset_download/WT_PPI_processing"

    logger.info("=" * 80)
    logger.info(f"RERUNNING {len(folders_with_duplicates)} FOLDERS WITH DUPLICATES")
    logger.info("=" * 80)

    idx_pdb_dict = parse_skempi()
    os.makedirs("logs", exist_ok=True)

    wt_folders = []
    base_folders = defaultdict(list)

    for edge_file_dir, row_index, results_path in folders_with_duplicates:
        rel_path = edge_file_dir.relative_to(results_path)
        parts = rel_path.parts

        if len(parts) >= 2 and parts[0] == "wt":
            pdb_code = parts[1]
            wt_folders.append((edge_file_dir, pdb_code))
        elif len(parts) >= 3:
            base_name = parts[0]
            folder_index = parts[2]
            base_folders[base_name].append((edge_file_dir, folder_index, row_index))

    if wt_folders:
        logger.info(f"Processing {len(wt_folders)} wildtype folders...")
        for edge_file_dir, pdb_code in tqdm(wt_folders, desc="Rerunning wildtype folders"):
            try:
                wt_input_path = Path(WT_PATH) / pdb_code
                if not wt_input_path.exists():
                    logger.error(f"Wildtype input path not found: {wt_input_path}")
                    continue

                halfs = sorted([d for d in os.listdir(wt_input_path) if d.startswith("H") and len(d) >= 4])
                if len(halfs) != 2:
                    logger.error(f"{pdb_code} is missing PDB halves.")
                    continue

                file1 = wt_input_path / halfs[0] / "hydrogensAdded.pdb"
                file2 = wt_input_path / halfs[1] / "hydrogensAdded.pdb"
                output = f"wt/{pdb_code}"

                logger.info(f"Rerunning wildtype: {pdb_code} -> {output}")
                calculate_diffbond(
                    str(file1),
                    str(file2),
                    output,
                    skempi_csv=None,
                    row_index=None,
                    skip_mutant=True,
                )
                logger.info(f"Successfully reran wildtype: {pdb_code}")
            except Exception as e:
                logger.error(f"Error rerunning wildtype folder {pdb_code}: {e}", exc_info=True)

    for base_name, folders in base_folders.items():
        try:
            base_num = int(base_name.split("-")[1])
        except (ValueError, IndexError):
            logger.error(f"Could not extract base number from {base_name}")
            continue

        base_input_path = Path(BASE_PATH_TEMPLATE.format(i=base_num))
        if not base_input_path.exists():
            logger.warning(f"Base input path does not exist: {base_input_path}, skipping")
            continue

        logger.info(f"Processing {len(folders)} folders in {base_name}...")
        for edge_file_dir, folder_index, row_index in tqdm(folders, desc=f"Rerunning {base_name}"):
            try:
                folder_path = base_input_path / folder_index
                if not folder_path.exists():
                    logger.error(f"Input folder not found: {folder_path}")
                    continue

                halfs = sorted([d for d in os.listdir(folder_path) if d.startswith("H") and len(d) >= 4])
                if len(halfs) != 2:
                    logger.error(f"{folder_index} is missing PDB halves.")
                    continue

                file1 = folder_path / halfs[0] / "half1.pdb"
                file2 = folder_path / halfs[1] / "half2.pdb"
                output = idx_pdb_dict.get(folder_index)
                if not output:
                    logger.error(f"Could not find output mapping for folder {folder_index}")
                    continue

                logger.info(f"Rerunning: {base_name}/{folder_index} -> {output}")
                calculate_diffbond(
                    str(file1),
                    str(file2),
                    output,
                    skempi_csv="datasets/MODIFIED-skempi_v2.csv" if not skip_mutant else None,
                    row_index=row_index if not skip_mutant else None,
                    skip_mutant=skip_mutant,
                )
                logger.info(f"Successfully reran: {base_name}/{folder_index}")
            except Exception as e:
                logger.error(f"Error rerunning folder {base_name}/{folder_index}: {e}", exc_info=True)

    logger.info("=" * 80)
    logger.info("FINISHED RERUNNING FOLDERS WITH DUPLICATES")
    logger.info("=" * 80)


def print_found_duplicates(
    results_root: Optional[Path] = None,
    test_folder: Optional[Path] = None,
) -> None:
    """Print detailed information about found duplicates without rerunning calculations."""
    logger.info("=" * 80)
    logger.info("DUPLICATE DETECTION: Finding and printing duplicate nodes")
    logger.info("=" * 80)

    folders_with_duplicates = find_folders_with_duplicates(results_root=results_root, test_folder=test_folder)

    if not folders_with_duplicates:
        logger.info("No folders with duplicates found.")
        return

    logger.info(f"\nFound {len(folders_with_duplicates)} folder(s) with duplicates:\n")

    for edge_file_dir, row_index, results_path in folders_with_duplicates:
        rel_path = edge_file_dir.relative_to(results_path) if results_path in edge_file_dir.parents else edge_file_dir
        logger.info(f"Folder: {rel_path}")
        if row_index:
            logger.info(f"  Row index: {row_index}")

        edge_files = list(edge_file_dir.glob("*_edges.txt"))
        edge_files = [
            f
            for f in edge_files
            if not (f.name == "mutant_edges.txt" and len(rel_path.parts) >= 3 and rel_path.parts[0] == "wt")
        ]

        found_duplicates_in_folder = False
        for edge_file in edge_files:
            duplicates = get_duplicates_from_file(edge_file)
            if duplicates:
                found_duplicates_in_folder = True
                logger.info(f"  File: {edge_file.name}")
                logger.info(f"    Found {len(duplicates)} duplicate position(s):")

                node_index, _, _ = read_edge_file(edge_file)
                for (chain, resseq), indices in sorted(duplicates.items()):
                    logger.info(f"      Position ({chain}, {resseq}):")
                    for idx in sorted(indices):
                        node = node_index.get(idx)
                        if node:
                            logger.info(f"        Index {idx}: {node}")

        if not found_duplicates_in_folder:
            logger.info("  (No duplicates found in edge files)")
        logger.info("")

    logger.info("=" * 80)
    logger.info(f"Total: {len(folders_with_duplicates)} folder(s) with duplicates")
    logger.info("=" * 80)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Detect duplicate nodes in edge files and rerun DiffBond calculations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--results-root",
        type=str,
        help="Root Results directory for duplicate detection (default: ./Results)",
    )

    parser.add_argument(
        "--skip-mutant",
        action="store_true",
        help="Skip mutant distance calculation even when SKEMPI CSV is available.",
    )

    parser.add_argument(
        "--print-only",
        action="store_true",
        help="Only print found duplicates without rerunning calculations.",
    )

    return parser.parse_args()


def main():
    """Main function for rerun-duplicates mode: detect duplicates and rerun calculations."""
    args = parse_arguments()

    os.makedirs("logs", exist_ok=True)

    results_root = Path(args.results_root) if args.results_root else None

    if args.print_only:
        print_found_duplicates(results_root=results_root, test_folder=None)
    else:
        logger.info("=" * 80)
        logger.info("RERUN DUPLICATES: Detecting folders with duplicate nodes")
        logger.info("=" * 80)

        logger.info("DiffBond instances will log to: logs/diffbond_all.log")
        logger.info("DiffBond failures will log to: logs/diffbond_failures.log")

        folders_with_duplicates = find_folders_with_duplicates(results_root=results_root, test_folder=None)

        if not folders_with_duplicates:
            logger.info("No folders with duplicates found. Nothing to rerun.")
            return

        logger.info(f"Found {len(folders_with_duplicates)} folders with duplicates")
        logger.info("Starting rerun process...")

        rerun_folders_with_duplicates(folders_with_duplicates=folders_with_duplicates, skip_mutant=args.skip_mutant)


if __name__ == "__main__":
    main()
