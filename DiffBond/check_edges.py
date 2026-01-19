#!/usr/bin/env python3
"""
Script to validate edge files and PDB files for correctness.

This script performs comprehensive validation:
1. Checks all PDB files with insertion codes have proper insertion codes added
2. Checks all mutant nodes have edges to all other nodes
3. Checks correct number of mutant_edges files in wt folder for each mutant in SKEMPI
4. Checks all edge files have proper headers: "Index,Node" and "Source,Target,Distance"
5. Checks each numbered folder has all required edge files
6. Checks no edge files are completely empty (should have headers at minimum)
7. Checks for duplicate nodes (same chain and normalized resSeq with different resName)
"""

import os
import csv
import sys
import re
import ast
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional, Any
from collections import defaultdict
from tqdm import tqdm

# Required edge file types that should exist in every numbered folder
REQUIRED_EDGE_FILES = [
    "cationpi_edges.txt",
    "contact_edges.txt",
    "covalent_edges.txt",
    "hbond_edges.txt",
    "ionic_edges.txt",
    "mutant_edges.txt",
    "saltbridge_hbond_edges.txt",
    "saltbridge_ionic_edges.txt",
]


def normalize_resseq(resseq: str) -> str:
    """Normalize residue sequence number by removing alternate location suffixes.

    Args:
        resseq: Residue sequence number string (may contain suffix like '100A')

    Returns:
        Normalized residue sequence number (numeric part only)
    """
    match = re.match(r"^(-?\d+)", str(resseq).strip())
    if match:
        return match.group(1)
    return str(resseq).strip()


def parse_node_string(node_str: str) -> Optional[Tuple[str, str, str]]:
    """Parse a node string representation to extract (chain, resSeq, resName).

    Args:
        node_str: String representation of node tuple, e.g., "('A', '60', 'ALA')"

    Returns:
        Tuple of (chain, resSeq, resName) or None if parsing fails
    """
    if not node_str or not node_str.strip():
        return None

    node_str = node_str.strip()
    # Remove surrounding quotes if present
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
    """Read an edge file and return nodes and edges.

    Args:
        edge_file: Path to edge file

    Returns:
        Tuple of (node_index_dict, edges_list, has_distance)
        - node_index_dict: Maps index -> (chain, resSeq, resName)
        - edges_list: List of [source, target] or [source, target, distance]
        - has_distance: Whether edges have distance column
    """
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

    # Find the edges header row
    edges_start_idx = None
    for i, row in enumerate(rows):
        if len(row) >= 2 and row[0].strip().lower() == "source" and row[1].strip().lower() == "target":
            edges_start_idx = i
            has_distance = len(row) >= 3 and row[2].strip().lower() == "distance"
            break

    if edges_start_idx is None:
        edges_start_idx = len(rows)

    # Parse nodes (format: Index,Node)
    # Skip header row if it exists
    start_row = 1 if rows[0][0].strip().lower() in ["index", "node"] else 0

    for i in range(start_row, edges_start_idx):
        if len(rows[i]) < 2:
            continue
        try:
            index = int(rows[i][0].strip())
            node_str = rows[i][1].strip()
            node_tuple = parse_node_string(node_str)
            if node_tuple:
                node_index[index] = node_tuple
        except (ValueError, IndexError):
            continue

    # Parse edges (format: Source,Target[,Distance])
    for i in range(edges_start_idx + 1, len(rows)):
        if len(rows[i]) < 2:
            continue
        try:
            source = int(rows[i][0].strip())
            target = int(rows[i][1].strip())
            if has_distance and len(rows[i]) >= 3:
                distance = float(rows[i][2].strip())
                edges.append([source, target, distance])
            else:
                edges.append([source, target])
        except (ValueError, IndexError):
            continue

    return node_index, edges, has_distance


def get_mutation_indices(edge_file_dir: Path) -> Set[int]:
    """Get set of mutation node indices from mutant_edges.txt.

    Args:
        edge_file_dir: Directory containing edge files

    Returns:
        Set of mutation node indices (Source indices from mutant_edges.txt)
    """
    mutation_indices = set()
    mutant_file = edge_file_dir / "mutant_edges.txt"

    if not mutant_file.exists():
        return mutation_indices

    try:
        node_index, edges, _ = read_edge_file(mutant_file)
        # Mutation nodes are the Source nodes in mutant edges
        for edge in edges:
            if len(edge) >= 1:
                mutation_indices.add(edge[0])  # Source index
    except Exception as e:
        print(f"Error reading mutation indices from {mutant_file}: {e}")

    return mutation_indices


def find_duplicate_nodes(node_index: Dict[int, Tuple[str, str, str]]) -> Dict[Tuple[str, str], List[int]]:
    """Find duplicate nodes where same (chain, resSeq) appears with different resName.

    This indicates that insertion codes were lost during initial calculation.
    For example: ('D', '30', 'THR') and ('D', '30', 'ASP') means '30' and '30a' were both
    written as '30', losing the insertion code.

    Args:
        node_index: Dictionary mapping index -> (chain, resSeq, resName)

    Returns:
        Dictionary mapping (chain, resSeq) -> list of indices with that position but different resName
    """
    position_to_indices = defaultdict(list)

    for index, (chain, resseq, resname) in node_index.items():
        # Use (chain, resSeq) as key - we want to find cases where same position
        # has different residue names (indicating lost insertion codes)
        position_key = (chain, resseq)
        position_to_indices[position_key].append((index, resname))

    # Filter to only positions with multiple nodes (potential duplicates)
    # But only flag as duplicates if they have different resName
    duplicates = {}
    for pos, index_resname_list in position_to_indices.items():
        if len(index_resname_list) > 1:
            # Check if there are different residue names at this position
            resnames = set(resname for _, resname in index_resname_list)
            if len(resnames) > 1:
                # Same (chain, resSeq) but different resName - this indicates lost insertion codes
                duplicates[pos] = [index for index, _ in index_resname_list]

    return duplicates


def detect_duplicates_in_file(
    edge_file: Path,
    mutation_indices: Set[int],
) -> bool:
    """Detect if an edge file has duplicate nodes (without fixing them).

    Args:
        edge_file: Path to edge file
        mutation_indices: Set of mutation node indices

    Returns:
        True if duplicates are found, False otherwise
    """
    node_index, edges, has_distance = read_edge_file(edge_file)

    if not node_index:
        return False

    duplicates = find_duplicate_nodes(node_index)
    return len(duplicates) > 0


def check_pdb_insertion_codes(pdb_file: Path) -> Tuple[bool, List[str]]:
    """Check if PDB file properly includes insertion codes in residue sequence numbers.

    Args:
        pdb_file: Path to PDB file

    Returns:
        Tuple of (is_valid, list_of_errors)
    """
    errors = []
    if not pdb_file.exists():
        return False, [f"PDB file does not exist: {pdb_file}"]

    with open(pdb_file, "r") as f:
        for line_num, line in enumerate(f, 1):
            if not line.startswith("ATOM") or len(line) < 27:
                continue

            insertion_code = line[26:27].strip()
            if insertion_code and not insertion_code.isalpha():
                errors.append(f"Line {line_num}: Invalid insertion code '{insertion_code}' (should be A-Z)")

            col26 = line[26]
            if col26 != " " and not col26.isalpha():
                errors.append(f"Line {line_num}: Column 26 has invalid character '{col26}'")

    return len(errors) == 0, errors


def check_mutant_edges_completeness(edge_file_dir: Path, results_path: Path) -> Tuple[bool, List[str]]:
    """Check that all mutant nodes have edges to all other nodes.

    Args:
        edge_file_dir: Directory containing edge files
        results_path: Root Results directory

    Returns:
        Tuple of (is_valid, list_of_errors)
    """
    # Skip this check for wt folders (ignore mutant_edges.txt in wt)
    is_wt_numbered = is_wt_numbered_folder(edge_file_dir, results_path)
    is_wt_pdb = is_wt_pdb_folder(edge_file_dir, results_path)
    if is_wt_numbered or is_wt_pdb:
        return True, []  # Skip check for wt folders

    errors = []
    mutant_file = edge_file_dir / "mutant_edges.txt"

    if not mutant_file.exists():
        return False, [f"mutant_edges.txt not found in {edge_file_dir}"]

    node_index, edges, _ = read_edge_file(mutant_file)

    if not node_index:
        return False, [f"mutant_edges.txt has no nodes in {edge_file_dir}"]

    # Get all mutant node indices (Source nodes in mutant edges)
    mutant_indices = {edge[0] for edge in edges if len(edge) >= 1}

    if not mutant_indices:
        return False, [f"No mutant nodes found in {edge_file_dir}/mutant_edges.txt"]

    # For each mutant node, check it has edges to all other nodes
    all_node_indices = set(node_index.keys())

    for mut_idx in mutant_indices:
        targets = {edge[1] for edge in edges if len(edge) >= 2 and edge[0] == mut_idx}
        expected_targets = all_node_indices - {mut_idx}
        missing_targets = expected_targets - targets

        if missing_targets:
            mut_node = node_index.get(mut_idx, "unknown")
            errors.append(
                f"Mutant node {mut_idx} ({mut_node}) missing edges to {len(missing_targets)} node(s): "
                f"{sorted(list(missing_targets))[:10]}{'...' if len(missing_targets) > 10 else ''}"
            )

    return len(errors) == 0, errors


def check_wt_mutant_edges_count(results_path: Path, skempi_csv: Path) -> Tuple[bool, List[str]]:
    """Check that wt folder has correct number of mutant_edges files for each PDB.

    Args:
        results_path: Root Results directory
        skempi_csv: Path to SKEMPI CSV file

    Returns:
        Tuple of (is_valid, list_of_errors)
    """
    errors = []
    wt_path = results_path / "wt"

    if not wt_path.exists():
        return True, []  # No wt folder, nothing to check

    if not skempi_csv.exists():
        return False, [f"SKEMPI CSV not found: {skempi_csv}"]

    # Read SKEMPI CSV to count mutations per PDB
    pdb_mutation_counts = defaultdict(int)
    with open(skempi_csv, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=";")
        for row in reader:
            if len(row) < 2 or not row[1].strip():
                continue
            pdb_code = row[0].strip().split("_")[0].upper()
            pdb_mutation_counts[pdb_code] += 1

    # Check each PDB folder in wt
    for pdb_dir in wt_path.iterdir():
        if not pdb_dir.is_dir():
            continue

        pdb_code = pdb_dir.name.upper()
        expected_count = pdb_mutation_counts.get(pdb_code, 0)

        if expected_count == 0:
            continue

        # Ignore mutant_edges.txt in pdb_dir itself (wt/PDB should not have it)
        # Only count mutant_edges.txt files in numbered subdirectories (wt/PDB/XXXXX)
        actual_count = sum(
            1
            for numbered_dir in pdb_dir.iterdir()
            if numbered_dir.is_dir() and numbered_dir.name.isdigit() and (numbered_dir / "mutant_edges.txt").exists()
        )

        if actual_count != expected_count:
            errors.append(f"PDB {pdb_code}: Expected {expected_count} mutant_edges.txt files, found {actual_count}")

    return len(errors) == 0, errors


def check_edge_file_headers(edge_file: Path) -> Tuple[bool, List[str]]:
    """Check that edge file has correct headers: "Index,Node" and "Source,Target,Distance".

    Args:
        edge_file: Path to edge file

    Returns:
        Tuple of (is_valid, list_of_errors)
    """
    errors = []

    if not edge_file.exists():
        return False, [f"Edge file does not exist: {edge_file}"]

    with open(edge_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        rows = list(reader)

    if not rows:
        return False, [f"Edge file is completely empty: {edge_file}"]

    # Check first header (should be "Index,Node")
    if len(rows[0]) < 2:
        return False, [f"First header row has insufficient columns: {edge_file}"]

    first_col = rows[0][0].strip().lower()
    second_col = rows[0][1].strip().lower()

    if not (first_col == "index" and second_col == "node"):
        errors.append(f"First header should be 'Index,Node' but found '{rows[0][0]},{rows[0][1]}' in {edge_file}")

    # Find edges header
    for i, row in enumerate(rows):
        if len(row) >= 2 and row[0].strip().lower() == "source" and row[1].strip().lower() == "target":
            has_distance = len(row) >= 3 and row[2].strip().lower() == "distance"

            # Check if there are actual edge rows
            if i + 1 < len(rows):
                has_edges = any(
                    len(rows[j]) >= 2 and rows[j][0].strip().isdigit() and rows[j][1].strip().isdigit()
                    for j in range(i + 1, len(rows))
                )
                if has_edges and not has_distance:
                    errors.append(f"Edges section has edges but missing 'Distance' column in {edge_file}")
            break

    return len(errors) == 0, errors


def is_wt_numbered_folder(edge_file_dir: Path, results_path: Path) -> bool:
    """Check if a folder is a numbered folder inside Results/wt/PDB/."""
    try:
        parts = edge_file_dir.relative_to(results_path).parts
        return len(parts) == 3 and parts[0] == "wt" and edge_file_dir.name.isdigit()
    except (ValueError, AttributeError):
        return False


def is_wt_pdb_folder(edge_file_dir: Path, results_path: Path) -> bool:
    """Check if a folder is the PDB folder itself in Results/wt/PDB."""
    try:
        parts = edge_file_dir.relative_to(results_path).parts
        return len(parts) == 2 and parts[0] == "wt"
    except (ValueError, AttributeError):
        return False


def check_required_edge_files(edge_file_dir: Path, results_path: Path) -> Tuple[bool, List[str]]:
    """Check that all required edge files exist in the directory.

    Special handling for wt folder structure:
    - wt/PDB/ folder: should have all files EXCEPT mutant_edges.txt
    - wt/PDB/XXXXX/ folder: should ONLY have mutant_edges.txt

    Args:
        edge_file_dir: Directory containing edge files
        results_path: Root Results directory

    Returns:
        Tuple of (is_valid, list_of_errors)
    """
    errors = []
    missing_files = []
    unexpected_files = []

    # Check if this is a numbered folder in wt/PDB structure
    is_wt_numbered = is_wt_numbered_folder(edge_file_dir, results_path)
    is_wt_pdb = is_wt_pdb_folder(edge_file_dir, results_path)

    if is_wt_numbered:
        # Numbered folder in wt/PDB should ONLY have mutant_edges.txt
        for required_file in REQUIRED_EDGE_FILES:
            file_path = edge_file_dir / required_file
            if required_file == "mutant_edges.txt":
                if not file_path.exists():
                    missing_files.append(required_file)
            else:
                if file_path.exists():
                    unexpected_files.append(required_file)

        if missing_files:
            errors.append(f"Missing mutant_edges.txt in {edge_file_dir}")
        if unexpected_files:
            errors.append(f"Unexpected edge files in wt numbered folder {edge_file_dir}: {', '.join(unexpected_files)}")

    elif is_wt_pdb:
        # PDB folder in wt should have all files EXCEPT mutant_edges.txt
        for required_file in REQUIRED_EDGE_FILES:
            if required_file == "mutant_edges.txt":
                continue  # Skip mutant_edges.txt for wt PDB folders
            file_path = edge_file_dir / required_file
            if not file_path.exists():
                missing_files.append(required_file)

        if missing_files:
            errors.append(f"Missing required edge files in wt PDB folder {edge_file_dir}: {', '.join(missing_files)}")

    else:
        # Regular folder (base-X/PDB/XXXXX) should have all required files
        for required_file in REQUIRED_EDGE_FILES:
            file_path = edge_file_dir / required_file
            if not file_path.exists():
                missing_files.append(required_file)

        if missing_files:
            errors.append(f"Missing required edge files in {edge_file_dir}: {', '.join(missing_files)}")

    return len(errors) == 0, errors


def check_empty_edge_files(edge_file_dir: Path, results_path: Path) -> Tuple[bool, List[str]]:
    """Check that no edge files are completely empty (should have headers at minimum).

    Args:
        edge_file_dir: Directory containing edge files
        results_path: Root Results directory

    Returns:
        Tuple of (is_valid, list_of_errors)
    """
    errors = []

    # Check which files should exist based on folder type
    is_wt_numbered = is_wt_numbered_folder(edge_file_dir, results_path)
    is_wt_pdb = is_wt_pdb_folder(edge_file_dir, results_path)

    files_to_check = []
    if is_wt_numbered:
        # Only check mutant_edges.txt
        files_to_check = ["mutant_edges.txt"]
    elif is_wt_pdb:
        # Check all files except mutant_edges.txt
        files_to_check = [f for f in REQUIRED_EDGE_FILES if f != "mutant_edges.txt"]
    else:
        # Check all required files
        files_to_check = REQUIRED_EDGE_FILES

    for required_file in files_to_check:
        file_path = edge_file_dir / required_file
        if not file_path.exists():
            continue  # Missing files are handled by check_required_edge_files

        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read().strip()

        if not content or not content.split("\n")[0].strip():
            errors.append(f"Edge file is completely empty or has no headers: {file_path}")

    return len(errors) == 0, errors


def validate_folder(
    edge_file_dir: Path,
    results_path: Path,
    skempi_csv: Optional[Path] = None,
    check_pdb: bool = False,
) -> Dict[str, Tuple[bool, List[str]]]:
    """Perform all validation checks on a folder.

    Args:
        edge_file_dir: Directory containing edge files
        results_path: Root Results directory
        skempi_csv: Path to SKEMPI CSV file (optional)
        check_pdb: Whether to check PDB files for insertion codes

    Returns:
        Dictionary mapping check_name -> (is_valid, list_of_errors)
    """
    results = {}

    # Check 1: Duplicate nodes
    mutation_indices = get_mutation_indices(edge_file_dir)
    has_duplicates = False
    duplicate_errors = []
    for edge_file_name in REQUIRED_EDGE_FILES:
        edge_file = edge_file_dir / edge_file_name
        if edge_file.exists():
            if detect_duplicates_in_file(edge_file, mutation_indices):
                has_duplicates = True
                duplicate_errors.append(f"Duplicate nodes found in {edge_file_name}")
    results["duplicates"] = (not has_duplicates, duplicate_errors)

    # Check 2: Mutant edges completeness
    results["mutant_edges_completeness"] = check_mutant_edges_completeness(edge_file_dir, results_path)

    # Check 3: Headers
    header_errors = []
    # Determine which files to check based on folder type
    is_wt_numbered = is_wt_numbered_folder(edge_file_dir, results_path)
    is_wt_pdb = is_wt_pdb_folder(edge_file_dir, results_path)

    files_to_check = []
    if is_wt_numbered:
        # Only check mutant_edges.txt
        files_to_check = ["mutant_edges.txt"]
    elif is_wt_pdb:
        # Check all files except mutant_edges.txt
        files_to_check = [f for f in REQUIRED_EDGE_FILES if f != "mutant_edges.txt"]
    else:
        # Check all required files
        files_to_check = REQUIRED_EDGE_FILES

    for edge_file_name in files_to_check:
        edge_file = edge_file_dir / edge_file_name
        if edge_file.exists():
            is_valid, errors = check_edge_file_headers(edge_file)
            if not is_valid:
                header_errors.extend(errors)
    results["headers"] = (len(header_errors) == 0, header_errors)

    # Check 4: Required files
    results["required_files"] = check_required_edge_files(edge_file_dir, results_path)

    # Check 5: Empty files
    results["empty_files"] = check_empty_edge_files(edge_file_dir, results_path)

    # Check 6: PDB insertion codes (if requested and PDB files exist)
    if check_pdb:
        pdb_errors = []
        pdb_dir = edge_file_dir / "pdb"
        if pdb_dir.exists():
            for pdb_file in pdb_dir.glob("*.pdb"):
                is_valid, errors = check_pdb_insertion_codes(pdb_file)
                if not is_valid:
                    pdb_errors.extend(errors)
        results["pdb_insertion_codes"] = (len(pdb_errors) == 0, pdb_errors)

    return results


def validate_all(
    results_root: Optional[Path] = None,
    skempi_csv: Optional[Path] = None,
    test_folder: Optional[Path] = None,
    check_pdb: bool = False,
) -> Dict[str, Any]:
    """Perform all validation checks across all folders.

    Args:
        results_root: Root Results directory
        skempi_csv: Path to SKEMPI CSV file
        test_folder: If provided, only check this specific folder
        check_pdb: Whether to check PDB files for insertion codes

    Returns:
        Dictionary with validation results
    """
    # Find Results directory
    if test_folder is not None:
        results_path = Path(test_folder)
        if not results_path.exists():
            print(f"ERROR: Test folder not found: {results_path}")
            return {}
    else:
        if results_root is None:
            results_path = Path("Results")
        else:
            results_path = Path(results_root)

        if not results_path.exists():
            print(f"ERROR: Results directory not found: {results_path}")
            return {}

    # Find SKEMPI CSV if not provided
    if skempi_csv is None:
        csv_path = Path("datasets/MODIFIED-skempi_v2.csv")
        if csv_path.exists():
            skempi_csv = csv_path

    # Find all folders to check
    if test_folder:
        # In test mode, check the folder itself
        folders_to_check = [results_path] if results_path.is_dir() else []
    else:
        # Find all folders to validate:
        # 1. Numbered folders in base-X/PDB/XXXXX
        # 2. Numbered folders in wt/PDB/XXXXX
        # 3. PDB folders in wt/PDB (to check they have all files except mutant_edges.txt)
        folders_to_check = []
        for base_dir in results_path.iterdir():
            if not base_dir.is_dir():
                continue

            for pdb_dir in base_dir.iterdir():
                if not pdb_dir.is_dir():
                    continue

                # Check if this is a wt PDB folder - add it to check list
                if base_dir.name == "wt":
                    # Check if it has any edge files (excluding mutant_edges.txt)
                    has_other_files = any(
                        (pdb_dir / f).exists() for f in REQUIRED_EDGE_FILES if f != "mutant_edges.txt"
                    )
                    if has_other_files:
                        folders_to_check.append(pdb_dir)

                # Check numbered subdirectories
                for numbered_dir in pdb_dir.iterdir():
                    if not numbered_dir.is_dir():
                        continue
                    # Check if it's a numbered directory (all digits)
                    if numbered_dir.name.isdigit():
                        # Check if it has any edge files
                        if any((numbered_dir / f).exists() for f in REQUIRED_EDGE_FILES):
                            folders_to_check.append(numbered_dir)

    print(f"Found {len(folders_to_check)} folder(s) to validate...")

    # Perform validation
    all_results = {
        "duplicates": [],
        "mutant_edges_completeness": [],
        "headers": [],
        "required_files": [],
        "empty_files": [],
    }

    if check_pdb:
        all_results["pdb_insertion_codes"] = []

    # Check wt mutant edges count (once for entire wt folder)
    if skempi_csv and results_path.exists():
        is_valid, errors = check_wt_mutant_edges_count(results_path, skempi_csv)
        all_results["wt_mutant_edges_count"] = (is_valid, errors)

    # Validate each folder
    for folder in tqdm(folders_to_check, desc="Validating folders"):
        results = validate_folder(folder, results_path, skempi_csv, check_pdb)

        for check_name, (is_valid, errors) in results.items():
            if not is_valid:
                all_results[check_name].append((folder, errors))

    return all_results


def main():
    """Main function to validate edge files and PDB files.

    Usage:
        python check_edges.py                          # Validate all edge files
        python check_edges.py --test <folder>         # Test on one folder only
        python check_edges.py --check-pdb             # Also check PDB insertion codes
    """
    # Parse arguments
    test_folder = None
    check_pdb = False

    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == "--test":
            if i + 1 >= len(sys.argv):
                print("Usage: python check_edges.py --test <folder_path>")
                print("Example: python check_edges.py --test Results/base-0/1ACB/00007")
                return
            test_folder = Path(sys.argv[i + 1])
            i += 2
        elif sys.argv[i] == "--check-pdb":
            check_pdb = True
            i += 1
        else:
            print(f"Unknown argument: {sys.argv[i]}")
            print("Usage: python check_edges.py [--test <folder>] [--check-pdb]")
            return

    print("=" * 60)
    if test_folder:
        print("TEST MODE: Validating single folder")
        print(f"Folder: {test_folder}")
    else:
        print("Validating all edge files")
    if check_pdb:
        print("(Including PDB insertion code checks)")
    print("=" * 60)

    # Get SKEMPI CSV path from environment or use default
    skempi_csv = os.environ.get("SKEMPI_CSV", None)
    if skempi_csv:
        skempi_csv = Path(skempi_csv)

    # Perform validation
    results = validate_all(
        results_root=None,
        skempi_csv=skempi_csv,
        test_folder=test_folder,
        check_pdb=check_pdb,
    )

    # Print summary
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)

    total_issues = 0

    # Check 1: Duplicates
    if "duplicates" in results:
        issues = results["duplicates"]
        if issues:
            total_issues += len(issues)
            print(f"\n❌ DUPLICATE NODES: {len(issues)} folder(s) with duplicates")
            for folder, errors in issues[:10]:
                print(f"  - {folder}")
                for error in errors[:3]:
                    print(f"    {error}")
            if len(issues) > 10:
                print(f"  ... and {len(issues) - 10} more")
        else:
            print("\n✅ DUPLICATE NODES: No duplicates found")

    # Check 2: Mutant edges completeness
    if "mutant_edges_completeness" in results:
        issues = results["mutant_edges_completeness"]
        if issues:
            total_issues += len(issues)
            print(f"\n❌ MUTANT EDGES COMPLETENESS: {len(issues)} folder(s) with incomplete mutant edges")
            for folder, errors in issues[:10]:
                print(f"  - {folder}")
                for error in errors[:3]:
                    print(f"    {error}")
            if len(issues) > 10:
                print(f"  ... and {len(issues) - 10} more")
        else:
            print("\n✅ MUTANT EDGES COMPLETENESS: All mutant nodes have edges to all other nodes")

    # Check 3: Headers
    if "headers" in results:
        issues = results["headers"]
        if issues:
            total_issues += len(issues)
            print(f"\n❌ HEADERS: {len(issues)} folder(s) with header issues")
            for folder, errors in issues[:10]:
                print(f"  - {folder}")
                for error in errors[:3]:
                    print(f"    {error}")
            if len(issues) > 10:
                print(f"  ... and {len(issues) - 10} more")
        else:
            print("\n✅ HEADERS: All edge files have correct headers")

    # Check 4: Required files
    if "required_files" in results:
        issues = results["required_files"]
        if issues:
            total_issues += len(issues)
            print(f"\n❌ REQUIRED FILES: {len(issues)} folder(s) missing required edge files")
            for folder, errors in issues[:10]:
                print(f"  - {folder}")
                for error in errors[:3]:
                    print(f"    {error}")
            if len(issues) > 10:
                print(f"  ... and {len(issues) - 10} more")
        else:
            print("\n✅ REQUIRED FILES: All folders have all required edge files")

    # Check 5: Empty files
    if "empty_files" in results:
        issues = results["empty_files"]
        if issues:
            total_issues += len(issues)
            print(f"\n❌ EMPTY FILES: {len(issues)} folder(s) with empty edge files")
            for folder, errors in issues[:10]:
                print(f"  - {folder}")
                for error in errors[:3]:
                    print(f"    {error}")
            if len(issues) > 10:
                print(f"  ... and {len(issues) - 10} more")
        else:
            print("\n✅ EMPTY FILES: No empty edge files found")

    # Check 6: PDB insertion codes
    if "pdb_insertion_codes" in results:
        issues = results["pdb_insertion_codes"]
        if issues:
            total_issues += len(issues)
            print(f"\n❌ PDB INSERTION CODES: {len(issues)} folder(s) with PDB insertion code issues")
            for folder, errors in issues[:10]:
                print(f"  - {folder}")
                for error in errors[:3]:
                    print(f"    {error}")
            if len(issues) > 10:
                print(f"  ... and {len(issues) - 10} more")
        else:
            print("\n✅ PDB INSERTION CODES: All PDB files have proper insertion codes")

    # Check 7: WT mutant edges count
    if "wt_mutant_edges_count" in results:
        is_valid, errors = results["wt_mutant_edges_count"]
        if not is_valid:
            total_issues += len(errors)
            print("\n❌ WT MUTANT EDGES COUNT: Issues found")
            for error in errors:
                print(f"  - {error}")
        else:
            print("\n✅ WT MUTANT EDGES COUNT: Correct number of mutant_edges files in wt folder")

    print("\n" + "=" * 60)
    if total_issues == 0:
        print("✅ ALL VALIDATION CHECKS PASSED")
    else:
        print(f"❌ TOTAL ISSUES FOUND: {total_issues}")
        print("\nTo rerun calculations for folders with duplicates, use:")
        print("  python rerun_duplicates.py")
    print("=" * 60)


if __name__ == "__main__":
    main()
