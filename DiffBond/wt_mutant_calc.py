"""
Calculate equivalent mutant edge measurements in wildtype structures.

For each mutation in the SKEMPI dataset, this script:
1. Identifies the equivalent wildtype residue (the original residue before mutation)
2. Calculates edges from that residue to all other nodes in the wildtype graph
3. Writes results to Results/wt/PDB/{SKEMPI_index}/mutant_edges.txt

This allows comparison of mutant edge patterns with their wildtype equivalents.

Command-line Parameters:
    --all
        Process all wildtype PDBs in the WT directory (automated mode)

    --pdb-code PDB
        PDB code (e.g., '1brs') - required if not using --all

    --pdb-file1 PATH
        Path to first PDB file (half1) - required if not using --all

    --pdb-file2 PATH
        Path to second PDB file (half2) - required if not using --all

    --wt-path PATH
        Path to wildtype directory (default: ../../SKEMPI_dataset_download/WT_PPI_processing)

    --skempi-csv PATH
        Path to SKEMPI CSV file (default: datasets/MODIFIED-skempi_v2.csv)

    -d, --distance FLOAT
        Distance threshold in Å

    --contact-distance {avg-atom,centroid}
        Distance calculation strategy (default: avg-atom)

    -v, --verbose
        Enable verbose logging

    --start-from PDB
        Start processing from this PDB code (inclusive, for --all mode)

    --end-at PDB
        End processing at this PDB code (exclusive, for --all mode)

    --test N
        Test mode: process only N PDBs (for --all mode)

    --start-mutant INDEX
        Start processing from this mutation index (inclusive)

    --end-mutant INDEX
        End processing at this mutation index (inclusive)

Example Commands:
    # Process all wildtype PDBs
    python wt_mutant_calc.py --all

    # Process a single PDB
    python wt_mutant_calc.py --pdb-code 1brs --pdb-file1 path/to/half1.pdb --pdb-file2 path/to/half2.pdb

    # Process with mutation range
    python wt_mutant_calc.py --all --start-mutant 5537 --end-mutant 5548

    # Test mode: process only 3 PDBs
    python wt_mutant_calc.py --all --test 3
"""

import os
import csv
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional
import re
from tqdm import tqdm

import src.utils.logging_handler as logging_handler
import src.utils.file_handler as file_handler
import src.utils.mutation_handler as mutation_handler
import src.utils.graph_handler as graph_handler
import src.core.interactions as interactions
from src.utils.residue_handler import build_residue_atoms
from src.utils.constants import MOLECULAR_MEASUREMENTS

logger = logging.getLogger(__name__)
failure_logger = logging.getLogger("failures")

AA_MAP = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}


def parse_mutation_token_for_wt(token: str) -> Optional[Tuple[str, str, str]]:
    """Parse a mutation token like 'AB32G' or 'RD100bA' to extract wildtype residue info.

    Returns:
        Tuple of (chain, resSeq, original_residue_name) or None if parsing fails
    """
    m = re.match(r"^([A-Z])([A-Z])(-?\d+[a-z]*)([A-Z])$", token)
    if not m:
        return None

    original_residue_code = m.group(1)
    chain = m.group(2)
    resseq = m.group(3)

    original_residue_name = AA_MAP.get(original_residue_code.upper())
    if not original_residue_name:
        return None

    return chain, resseq, original_residue_name


def get_all_mutations_for_pdb(
    skempi_csv: Path,
    pdb_code: str,
    start_mutant: Optional[int] = None,
    end_mutant: Optional[int] = None,
) -> Dict[int, List[str]]:
    """Get all mutations for a given PDB code from SKEMPI CSV."""
    mutations_by_index: Dict[int, List[str]] = {}

    try:
        with open(skempi_csv, "r") as f:
            reader = csv.reader(f, delimiter=";")
            for i, row in enumerate(reader):
                row_index = i + 1

                if start_mutant is not None and row_index < start_mutant:
                    continue
                if end_mutant is not None and row_index > end_mutant:
                    continue

                if not row:
                    continue

                pdb_from_row = row[0].split("_")[0] if "_" in row[0] else row[0]
                if pdb_from_row.upper() != pdb_code.upper():
                    continue

                if len(row) > 1 and row[1].strip():
                    tokens = [t.strip() for t in row[1].strip().split(",") if t.strip()]
                    if tokens:
                        mutations_by_index[row_index] = tokens
    except Exception as err:
        failure_logger.exception(f"[{pdb_code}] Failed reading SKEMPI CSV: {err}")

    return mutations_by_index


def compute_wt_equivalent_edges(
    pdb_data: List[List[List[str]]],
    mut_specs: List[str],
    union_residue_keys: Set[Tuple[str, str, str]],
    contact_distance_strategy: str,
    pdb_code: str,
) -> Tuple[Dict[int, Tuple[str, str, str]], List[Tuple[int, int, float]]]:
    """Compute distances from wildtype equivalent residues to all other residue nodes."""
    res_atoms_all: Dict[Tuple[str, str, str], List[List[str]]] = {}
    for entity in pdb_data:
        res_atoms_all.update(build_residue_atoms(entity))

    chain_res_to_key: Dict[Tuple[str, str], Tuple[str, str, str]] = {}
    for key in res_atoms_all.keys():
        chain_res_to_key[(str(key[0]).strip(), str(key[1]).strip())] = key

    wt_equiv_keys_list: List[Tuple[str, str, str]] = []
    for token in mut_specs:
        parsed = parse_mutation_token_for_wt(token)
        if not parsed:
            failure_logger.warning(f"[{pdb_code}] Failed to parse mutation token '{token}'. Skipping.")
            continue

        chain, resseq, expected_resname = parsed

        mut_key = None
        for (ch, res), key in chain_res_to_key.items():
            if ch.upper() == chain.upper() and res.upper() == resseq.upper():
                mut_key = key
                break

        if mut_key is None:
            available_chains = (
                ", ".join(sorted(set(ch for ch, _ in chain_res_to_key.keys()))) if chain_res_to_key else "none"
            )
            failure_logger.warning(
                f"[{pdb_code}] Could not find wildtype equivalent residue for token '{token}' "
                f"(chain={chain}, resSeq={resseq}). Available chains: {available_chains}. Skipping."
            )
            continue

        actual_resname = mut_key[2].strip().upper()
        expected_resname_upper = expected_resname.upper()

        if actual_resname != expected_resname_upper:
            failure_logger.warning(
                f"[{pdb_code}] Wildtype residue name mismatch for token '{token}': "
                f"expected {expected_resname_upper}, found {actual_resname} at {mut_key}. Using found residue."
            )

        wt_atoms = res_atoms_all.get(mut_key)
        if not wt_atoms:
            failure_logger.warning(
                f"[{pdb_code}] No atoms found for wildtype equivalent residue {mut_key} from token '{token}'. Skipping."
            )
            continue

        wt_equiv_keys_list.append(mut_key)

    if not wt_equiv_keys_list:
        failure_logger.warning(f"[{pdb_code}] No valid wildtype equivalent residues found for mutations: {mut_specs}")
        return {}, []

    wt_index_map: Dict[int, Tuple[str, str, str]] = {}
    wt_inv_map: Dict[Tuple[str, str, str], int] = {}
    current_idx = 0

    for wt_key in wt_equiv_keys_list:
        if wt_key not in wt_inv_map:
            wt_index_map[current_idx] = wt_key
            wt_inv_map[wt_key] = current_idx
            current_idx += 1

    other_residue_keys = sorted(union_residue_keys - set(wt_equiv_keys_list))
    for other_key in other_residue_keys:
        if other_key not in wt_inv_map:
            wt_index_map[current_idx] = other_key
            wt_inv_map[other_key] = current_idx
            current_idx += 1

    all_wt_edges: List[Tuple[int, int, float]] = []
    for wt_key in wt_equiv_keys_list:
        wt_atoms = res_atoms_all.get(wt_key)
        if not wt_atoms:
            continue

        wt_idx = wt_inv_map.get(wt_key)
        if wt_idx is None:
            continue

        all_other_keys = set(wt_index_map.values())
        all_other_keys.discard(wt_key)

        for other_key in sorted(all_other_keys):
            other_atoms = res_atoms_all.get(other_key)
            if not other_atoms:
                continue

            other_idx = wt_inv_map.get(other_key)
            if other_idx is None:
                continue

            dval = mutation_handler.residue_distance(wt_atoms, other_atoms, contact_distance_strategy)
            all_wt_edges.append((wt_idx, other_idx, dval))

    return wt_index_map, all_wt_edges


def process_wt_pdb_with_mutations(
    pdb_code: str,
    pdb_file1: Path,
    pdb_file2: Path,
    skempi_csv: Path,
    distance: float = MOLECULAR_MEASUREMENTS["default_max_distance"],
    contact_distance_strategy: str = "avg-atom",
    verbose: bool = False,
    start_mutant: Optional[int] = None,
    end_mutant: Optional[int] = None,
):
    """Process a wildtype PDB and calculate equivalent edges for all its mutations."""
    try:
        pdb_data1 = file_handler.parse_PDB_file(pdb_file1)
        pdb_data2 = file_handler.parse_PDB_file(pdb_file2)
        if not pdb_data1 or not pdb_data2:
            failure_logger.error(f"[{pdb_code}] No data found in PDB files")
            return
        pdb_data = [pdb_data1, pdb_data2]
    except Exception as e:
        failure_logger.exception(f"[{pdb_code}] Failed to parse PDB files: {e}")
        return

    mutations_by_index = get_all_mutations_for_pdb(skempi_csv, pdb_code, start_mutant, end_mutant)
    if not mutations_by_index:
        return

    logger.info(f"Calculating wildtype interface contact residues for {pdb_code}...")
    try:
        _, _, res_residue_edges = interactions.c_mode(
            pdb_data,
            distance,
            verbose=False,
            node_level="residue",
            contact_distance_strategy=contact_distance_strategy,
        )
    except Exception as e:
        failure_logger.exception(f"[{pdb_code}] Failed to calculate contact residues: {e}")
        return

    union_residue_keys: Set[Tuple[str, str, str]] = set()
    for r1, r2, _d in res_residue_edges:
        union_residue_keys.add(r1)
        union_residue_keys.add(r2)

    mutations_by_wt_residue: Dict[Tuple[str, str], List[Tuple[int, List[str]]]] = {}

    for row_index, mut_specs in mutations_by_index.items():
        for token in mut_specs:
            parsed = parse_mutation_token_for_wt(token)
            if parsed:
                chain, resseq, _ = parsed
                wt_residue_key = (chain, resseq)
                if wt_residue_key not in mutations_by_wt_residue:
                    mutations_by_wt_residue[wt_residue_key] = []
                mutations_by_wt_residue[wt_residue_key].append((row_index, mut_specs))

    logger.info(
        f"Grouped {len(mutations_by_index)} mutation entries into {len(mutations_by_wt_residue)} "
        f"unique wildtype residues (caching enabled)"
    )

    edge_cache: Dict[Tuple[str, str], Tuple[Dict[int, Tuple[str, str, str]], List[Tuple[int, int, float]]]] = {}

    for wt_residue_key, mutation_entries in mutations_by_wt_residue.items():
        chain, resseq = wt_residue_key

        _, mut_specs_all = mutation_entries[0]
        mut_specs = [
            token
            for token in mut_specs_all
            if parse_mutation_token_for_wt(token) and parse_mutation_token_for_wt(token)[:2] == wt_residue_key
        ]

        if not mut_specs:
            failure_logger.warning(
                f"[{pdb_code}] No valid mutations found for residue {chain}{resseq} "
                f"in mutation entry {mutation_entries[0][0]}"
            )
            continue

        try:
            wt_index_map, wt_edges = compute_wt_equivalent_edges(
                pdb_data,
                mut_specs,
                union_residue_keys,
                contact_distance_strategy,
                pdb_code,
            )

            if wt_edges:
                edge_cache[wt_residue_key] = (wt_index_map, wt_edges)
                logger.info(
                    f"Computed and cached edge graph for {chain}{resseq}: "
                    f"{len(wt_edges)} edges from {len(set(src for src, _, _ in wt_edges))} node(s)"
                )
            else:
                failure_logger.warning(
                    f"[{pdb_code}] No wildtype equivalent edges computed for "
                    f"residue {chain}{resseq}. Mutation token(s): {mut_specs}."
                )
        except Exception as err:
            failure_logger.exception(
                f"[{pdb_code}] Wildtype equivalent distance computation failed for " f"residue {chain}{resseq}: {err}"
            )
            continue

    for wt_residue_key, mutation_entries in mutations_by_wt_residue.items():
        if wt_residue_key not in edge_cache:
            continue

        wt_index_map, wt_edges = edge_cache[wt_residue_key]
        chain, resseq = wt_residue_key

        for row_index, mut_specs in mutation_entries:
            formatted_index = f"{row_index:05d}"
            output_dir = Path("Results") / "wt" / pdb_code / formatted_index
            output_dir.mkdir(parents=True, exist_ok=True)

            try:
                graph_handler.write_index_edge_mapping(
                    index=wt_index_map,
                    edges=wt_edges,
                    name="mutant",
                    output_dir=output_dir,
                )
                num_equiv = len(set(src for src, _, _ in wt_edges))
                logger.info(
                    f"Wrote wildtype equivalent edge graph for {pdb_code}/{formatted_index} "
                    f"(residue {chain}{resseq}, cached): "
                    f"{len(wt_edges)} edges from {num_equiv} wildtype equivalent node(s)"
                )
            except Exception as e:
                failure_logger.exception(
                    f"[{pdb_code}] Failed writing wildtype equivalent edge graph for "
                    f"{pdb_code}/{formatted_index}: {e}"
                )


def find_wt_pdb_folders(wt_path: Path) -> Dict[str, Tuple[Path, Path]]:
    """Find all wildtype PDB folders and their half files."""
    pdb_folders: Dict[str, Tuple[Path, Path]] = {}

    if not wt_path.exists():
        logger.error(f"Wildtype path does not exist: {wt_path}")
        return pdb_folders

    for pdb_folder in wt_path.iterdir():
        if not pdb_folder.is_dir():
            continue

        pdb_code = pdb_folder.name
        halfs = sorted([d for d in pdb_folder.iterdir() if d.is_dir() and d.name.startswith("H") and len(d.name) >= 4])

        if len(halfs) != 2:
            logger.warning(f"{pdb_code} is missing PDB halves (found {len(halfs)} halves)")
            continue

        file1 = halfs[0] / "hydrogensAdded.pdb"
        file2 = halfs[1] / "hydrogensAdded.pdb"

        if not file1.exists() or not file2.exists():
            logger.warning(f"{pdb_code} is missing hydrogensAdded.pdb files")
            continue

        pdb_folders[pdb_code] = (file1, file2)

    return pdb_folders


def process_all_wt_pdbs(
    wt_path: Path,
    skempi_csv: Path,
    distance: float = MOLECULAR_MEASUREMENTS["default_max_distance"],
    contact_distance_strategy: str = "avg-atom",
    verbose: bool = False,
    start_from: Optional[str] = None,
    end_at: Optional[str] = None,
    test_samples: Optional[int] = None,
    start_mutant: Optional[int] = None,
    end_mutant: Optional[int] = None,
):
    """Process all wildtype PDBs and calculate equivalent edges for their mutations."""
    logger.info("=" * 60)
    logger.info("Starting automated wildtype equivalent edge calculation")
    logger.info(f"Wildtype path: {wt_path}")
    logger.info(f"SKEMPI CSV: {skempi_csv}")
    logger.info("=" * 60)

    pdb_folders = find_wt_pdb_folders(wt_path)
    if not pdb_folders:
        failure_logger.error("No valid wildtype PDB folders found")
        return

    pdb_codes = sorted(pdb_folders.keys())

    if start_from:
        pdb_codes = [pdb for pdb in pdb_codes if pdb >= start_from]

    if end_at:
        pdb_codes = [pdb for pdb in pdb_codes if pdb < end_at]

    if test_samples:
        original_count = len(pdb_codes)
        pdb_codes = pdb_codes[:test_samples]
        logger.info(f"TEST MODE: Processing {len(pdb_codes)} out of {original_count} PDBs")

    total_pdbs = len(pdb_codes)
    logger.info(f"Found {total_pdbs} PDBs to process")

    with tqdm(total=total_pdbs, desc="Processing wildtype PDBs") as progress_bar:
        for pdb_code in pdb_codes:
            try:
                file1, file2 = pdb_folders[pdb_code]
                process_wt_pdb_with_mutations(
                    pdb_code=pdb_code,
                    pdb_file1=file1,
                    pdb_file2=file2,
                    skempi_csv=skempi_csv,
                    distance=distance,
                    contact_distance_strategy=contact_distance_strategy,
                    verbose=verbose,
                    start_mutant=start_mutant,
                    end_mutant=end_mutant,
                )
            except Exception as exc:
                failure_logger.exception(f"[{pdb_code}] Error processing {pdb_code}: {exc}")
            finally:
                progress_bar.update(1)

    logger.info("=" * 60)
    logger.info("Finished processing all wildtype PDBs")
    logger.info("=" * 60)


def main():
    """Main entry point for processing wildtype PDBs with mutation equivalents."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Calculate equivalent mutant edge measurements in wildtype structures.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("--all", action="store_true", help="Process all wildtype PDBs in the WT directory")

    parser.add_argument("--pdb-code", type=str, help="PDB code (e.g., '1brs') - required if not using --all")

    parser.add_argument("--pdb-file1", type=Path, help="Path to first PDB file (half1) - required if not using --all")

    parser.add_argument("--pdb-file2", type=Path, help="Path to second PDB file (half2) - required if not using --all")

    parser.add_argument(
        "--wt-path",
        type=Path,
        default=Path("../../SKEMPI_dataset_download/WT_PPI_processing"),
        help="Path to wildtype directory",
    )

    parser.add_argument(
        "--skempi-csv",
        type=Path,
        default=Path("datasets/MODIFIED-skempi_v2.csv"),
        help="Path to SKEMPI CSV file",
    )

    parser.add_argument(
        "-d",
        "--distance",
        type=float,
        default=MOLECULAR_MEASUREMENTS["default_max_distance"],
        help="Distance threshold in Å",
    )

    parser.add_argument(
        "--contact-distance",
        choices=["avg-atom", "centroid"],
        default="avg-atom",
        help="Distance calculation strategy",
    )

    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")

    parser.add_argument(
        "--start-from",
        type=str,
        help="Start processing from this PDB code (inclusive, for --all mode)",
    )

    parser.add_argument("--end-at", type=str, help="End processing at this PDB code (exclusive, for --all mode)")

    parser.add_argument("--test", type=int, metavar="N", help="Test mode: process only N PDBs (for --all mode)")

    parser.add_argument(
        "--start-mutant",
        type=int,
        metavar="INDEX",
        help="Start processing from this mutation index (inclusive)",
    )

    parser.add_argument(
        "--end-mutant",
        type=int,
        metavar="INDEX",
        help="End processing at this mutation index (inclusive)",
    )

    args = parser.parse_args()

    os.makedirs("logs", exist_ok=True)
    logging_handler.setup_logging(verbose=args.verbose)

    if args.all:
        if not args.wt_path.exists():
            failure_logger.error(f"Wildtype path does not exist: {args.wt_path}")
            return

        if not args.skempi_csv.exists():
            failure_logger.error(f"SKEMPI CSV file does not exist: {args.skempi_csv}")
            return

        process_all_wt_pdbs(
            wt_path=args.wt_path,
            skempi_csv=args.skempi_csv,
            distance=args.distance,
            contact_distance_strategy=args.contact_distance,
            verbose=args.verbose,
            start_from=args.start_from,
            end_at=args.end_at,
            test_samples=args.test,
            start_mutant=args.start_mutant,
            end_mutant=args.end_mutant,
        )
    else:
        if not args.pdb_code or not args.pdb_file1 or not args.pdb_file2:
            parser.error("Either --all must be specified, or --pdb-code, --pdb-file1, and --pdb-file2 must be provided")

        if not args.skempi_csv.exists():
            failure_logger.error(f"SKEMPI CSV file does not exist: {args.skempi_csv}")
            return

        process_wt_pdb_with_mutations(
            pdb_code=args.pdb_code,
            pdb_file1=args.pdb_file1,
            pdb_file2=args.pdb_file2,
            skempi_csv=args.skempi_csv,
            distance=args.distance,
            contact_distance_strategy=args.contact_distance,
            verbose=args.verbose,
            start_mutant=args.start_mutant,
            end_mutant=args.end_mutant,
        )


if __name__ == "__main__":
    main()
