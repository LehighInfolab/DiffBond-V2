"""
Utility functions for processing SKEMPI mutations and computing mutant-to-residue distances.
"""

import csv
import re
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Set, Any, Optional
from math import sqrt

from src.utils.residue_handler import normalize_resseq, residue_key, build_residue_atoms
from src.utils.distance import centroid_of, residue_distance
from src.utils.constants import SINGLE_TO_THREE_LETTER

logger = logging.getLogger(__name__)
failure_logger = logging.getLogger("failures")


def parse_mutation_token(token: str) -> Optional[Tuple[str, str, str, str]]:
    """Parse a mutation token like 'LI45G' or 'RD100bA' into (original_wt, chain, resSeq, mutated).

    Args:
        token: Mutation token in format like 'LI45G' or 'RD100bA'
               Format: original_residue + chain + residue_number + mutated_residue
               Residue number can be digits (e.g., '45'), digits with letter (e.g., '100b', '60a'),
               or negative (e.g., '-25')

    Returns:
        Tuple of (original_wt, chain, resSeq, mutated) or None if parsing fails
    """
    # Match: original_residue (1 letter) + chain (1 letter) + residue_number
    # (optional negative, digits, optional lowercase letters) + mutated_residue (1 letter)
    m = re.match(r"^([A-Z])([A-Z])(-?\d+[a-z]*)([A-Z])$", token)
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3), m.group(4)


def read_skempi_mutations(skempi_csv: Path, row_index: int) -> List[str]:
    """Read mutation specifications from a SKEMPI CSV file.

    CSV format: row[0] = PDB code, row[1] = mutation tokens (correct position), row[2] = literature position (not used)

    Args:
        skempi_csv: Path to SKEMPI CSV file (semicolon-delimited, no header row)
        row_index: Row index from folder name (e.g., folder "00001" -> row_index=1)

    Returns:
        List of mutation token strings (may be comma-separated in CSV, will be split)
        Note: Uses row[1] which contains the correct mutation position, not row[2] (literature position)
    """
    target_row = row_index - 1
    mut_specs: List[str] = []
    with open(skempi_csv, "r") as f:
        reader = csv.reader(f, delimiter=";")
        for i, row in enumerate(reader):
            if i == target_row:
                # Use row[1] for mutations (correct position), not row[2] (literature position only)
                if len(row) > 1 and row[1].strip():
                    tokens = [t.strip() for t in row[1].strip().split(",") if t.strip()]
                    mut_specs.extend(tokens)
                break
    return mut_specs


def _is_residue_edge_format(edge_item: Any) -> bool:
    """Check if an edge item is in residue-level format.

    Residue-level format: (Tuple[str, str, str], Tuple[str, str, str], float)
    Atom-level format: (List[str], List[str])

    Args:
        edge_item: An edge from raw_edges list

    Returns:
        True if residue-level format (length 3 and first element is a 3-tuple), False otherwise
    """
    if not isinstance(edge_item, (list, tuple)) or len(edge_item) != 3:
        return False
    if not isinstance(edge_item[0], (list, tuple)) or len(edge_item[0]) != 3:
        return False
    return all(isinstance(x, str) for x in edge_item[0])


def _extract_residue_keys_from_edges(raw_edges: List[Any]) -> Set[Tuple[str, str, str]]:
    """Extract residue keys from raw edges, handling both atom-level and residue-level formats.

    Args:
        raw_edges: List of raw edges (either atom pairs or residue triples)

    Returns:
        Set of residue keys extracted from the edges
    """
    keys: Set[Tuple[str, str, str]] = set()

    if not raw_edges:
        return keys

    # Check format of first edge to determine structure
    is_residue_format = _is_residue_edge_format(raw_edges[0])

    if is_residue_format:
        # Residue-level: (resKey, resKey, distance)
        for edge in raw_edges:
            if isinstance(edge, (list, tuple)) and len(edge) >= 2:
                r1, r2 = edge[0], edge[1]
                if isinstance(r1, (list, tuple)) and len(r1) == 3:
                    keys.add(tuple(r1))
                if isinstance(r2, (list, tuple)) and len(r2) == 3:
                    keys.add(tuple(r2))
    else:
        # Atom-level: (atom1, atom2) where atoms are lists
        for edge in raw_edges:
            if isinstance(edge, (list, tuple)) and len(edge) >= 2:
                a1, a2 = edge[0], edge[1]
                if isinstance(a1, list):
                    keys.add(residue_key(a1))
                if isinstance(a2, list):
                    keys.add(residue_key(a2))

    return keys


def collect_residue_keys_from_edge_dict(
    edge_dict: Dict[str, List[Any]],
    res_residue_edges: List[Tuple],
) -> Set[Tuple[str, str, str]]:
    """Collect all residue keys appearing in interaction outputs.

    Args:
        edge_dict: Dictionary mapping interaction type to [index, edges, raw_edges]
        res_residue_edges: Residue-level contact edges

    Returns:
        Set of residue keys (chain, resSeq, resName)
    """
    union_residue_keys: Set[Tuple[str, str, str]] = set()

    # From contact residue edges
    for r1, r2, _d in res_residue_edges:
        union_residue_keys.add(r1)
        union_residue_keys.add(r2)

    # From ionic (if present)
    if "ionic" in edge_dict and edge_dict["ionic"][2]:
        union_residue_keys.update(_extract_residue_keys_from_edges(edge_dict["ionic"][2]))

    # From hbond (if present)
    if "hbond" in edge_dict and edge_dict["hbond"][2]:
        union_residue_keys.update(_extract_residue_keys_from_edges(edge_dict["hbond"][2]))

    # From salt bridges
    if "saltbridge_ionic" in edge_dict or "saltbridge_hbond" in edge_dict:
        for sb_type in ["saltbridge_ionic", "saltbridge_hbond"]:
            if sb_type in edge_dict:
                sb_raw = edge_dict[sb_type][2]
                if sb_raw:
                    union_residue_keys.update(_extract_residue_keys_from_edges(sb_raw))

    # From adjacent (future mode), if present
    if "adjacent" in edge_dict and edge_dict["adjacent"][2]:
        union_residue_keys.update(_extract_residue_keys_from_edges(edge_dict["adjacent"][2]))

    return union_residue_keys


def compute_mutant_distances(
    pdb_data: List[List[List[str]]],
    mut_specs: List[str],
    union_residue_keys: Set[Tuple[str, str, str]],
    contact_distance_strategy: str,
) -> Tuple[Dict[int, Tuple[str, str, str]], List[Tuple[int, int, float]]]:
    """Compute distances from mutant residues to all other residue nodes.

    Each mutation gets a unique node index (0, 1, 2, ...), and each mutation node
    connects to all other nodes (other mutations and residue nodes from interactions).

    Args:
        pdb_data: List of PDB data entities
        mut_specs: List of mutation token strings
        union_residue_keys: Set of all residue keys to compute distances to
        contact_distance_strategy: Distance calculation strategy ('avg-atom' or 'centroid')

    Returns:
        Tuple of (index_map, edges) where index_map maps indices to residue keys
        and edges are (mutant_idx, other_idx, distance) tuples
    """
    logger.info(f"Starting mutation edge calculation for {len(mut_specs)} mutation token(s): {mut_specs}")

    # Build residue-to-atoms mapping for all entities
    res_atoms_all: Dict[Tuple[str, str, str], List[List[str]]] = {}
    for entity in pdb_data:
        res_atoms_all.update(build_residue_atoms(entity))

    # First pass: parse and validate all mutations, collect their keys
    mut_keys_list: List[Tuple[str, str, str]] = []
    for token in mut_specs:
        parsed = parse_mutation_token(token)
        if not parsed:
            failure_logger.warning(
                f"Failed to parse mutation token '{token}'. "
                f"Expected format: single letter + chain + residue number + single letter (e.g., 'LI45G'). "
                f"Skipping this mutation."
            )
            continue

        original_wt, chain, resseq, mutated = parsed
        # Convert single-letter code to three-letter code for comparison
        mutated_three_letter = SINGLE_TO_THREE_LETTER.get(mutated.upper())
        if not mutated_three_letter:
            failure_logger.warning(
                f"Unknown amino acid code '{mutated}' in mutation token '{token}'. Skipping this mutation."
            )
            continue

        # Check if resseq has a letter suffix (e.g., "60b", "100a")
        has_suffix = normalize_resseq(resseq) != resseq
        normalized_resseq = normalize_resseq(resseq)

        mut_key = None

        # Priority 1: If resseq has a suffix (e.g., "100a"), try to find that exact variant first
        if has_suffix:
            for key in res_atoms_all.keys():
                ch, res, resname = str(key[0]).strip(), str(key[1]).strip(), key[2]
                if (
                    ch.upper() == chain.upper()
                    and res.upper() == resseq.upper()  # Exact match including suffix
                    and resname == mutated_three_letter
                ):
                    mut_key = key
                    logger.debug(
                        f"Found exact variant {res} for mutation token '{token}': "
                        f"chain={ch}, mutated residue={mutated} ({mutated_three_letter})"
                    )
                    break

        # Priority 2: Try ANY variant at the normalized position that matches the mutant residue name
        # This handles cases where mutation specifies "100a" but only "100b" has the mutant residue
        if mut_key is None:
            for key in res_atoms_all.keys():
                ch, res, resname = str(key[0]).strip(), str(key[1]).strip(), key[2]
                if (
                    ch.upper() == chain.upper()
                    and normalize_resseq(res) == normalized_resseq
                    and resname == mutated_three_letter  # Match mutant residue name (most important)
                ):
                    # CRITICAL: Use the actual key from res_atoms_all to preserve insertion codes
                    # Do NOT create a new key based on the mutation token's resseq
                    mut_key = key
                    logger.info(
                        f"Found variant {res} at position {normalized_resseq} (mutation specified {resseq}) "
                        f"for mutation token '{token}': chain={ch}, mutated residue={mutated} ({mutated_three_letter}), "
                        f"using key {mut_key} from res_atoms_all"
                    )
                    break

        # If still not found, provide helpful diagnostics
        if mut_key is None:
            # Check if residue number exists on any chain
            matching_residues = []
            for key in res_atoms_all.keys():
                ch, res = str(key[0]).strip(), str(key[1]).strip()
                if normalize_resseq(res) == normalized_resseq:
                    matching_residues.append((ch, res, key))

            # Get all unique chains in the PDB for reference
            all_chains = sorted(set(str(key[0]).strip() for key in res_atoms_all.keys()))

            if matching_residues:
                # Residue number exists but on different chain(s) or wrong variant
                available_chains = sorted(set(ch for ch, _, _ in matching_residues))
                available_variants = sorted(set(res for _, res, _ in matching_residues))
                failure_logger.warning(
                    f"Could not find mutant residue for token '{token}' (chain={chain}, resSeq={resseq}) "
                    f"in PDB structure. Found variants {', '.join(available_variants)} at position {normalized_resseq} "
                    f"on chain(s): {', '.join(available_chains)}. "
                    f"Available chains in PDB: {', '.join(all_chains) if all_chains else 'none'}. "
                    f"The chain identifier in the mutation token may be incorrect or the residue "
                    f"may be in a different chain in the half PDB files."
                )
            else:
                # Residue number doesn't exist at all
                failure_logger.warning(
                    f"Could not find mutant residue for token '{token}' (chain={chain}, resSeq={resseq}) "
                    f"in PDB structure. Available chains in PDB: "
                    f"{', '.join(all_chains) if all_chains else 'none'}. "
                    f"The residue may not exist in the input PDB files "
                    f"(half1.pdb/half2.pdb may have different numbering)."
                )
            continue

        mut_atoms = res_atoms_all.get(mut_key)
        if not mut_atoms:
            failure_logger.warning(
                f"No atoms found for mutant residue {mut_key} from token '{token}'. "
                f"This may indicate a data parsing issue."
            )
            continue

        mut_keys_list.append(mut_key)

    # Warn if no mutations were successfully found due to missing variants
    if not mut_keys_list:
        failure_logger.warning(
            f"Mutation edge calculation cannot be run: no mutations were found in PDB structure. "
            f"Mutation token(s) provided: {mut_specs}. "
            f"All mutations were skipped due to missing variants or invalid tokens."
        )
        return {}, []

    # Build index mapping: mutations get indices 0, 1, 2, ..., then other residues continue
    mut_index_map: Dict[int, Tuple[str, str, str]] = {}
    mut_inv_map: Dict[Tuple[str, str, str], int] = {}
    current_idx = 0

    # Index all mutations first (0, 1, 2, ...)
    for mut_key in mut_keys_list:
        if mut_key not in mut_inv_map:
            mut_index_map[current_idx] = mut_key
            mut_inv_map[mut_key] = current_idx
            current_idx += 1

    # Index all other residue nodes (excluding mutations themselves)
    other_residue_keys = sorted(union_residue_keys - set(mut_keys_list))
    for other_key in other_residue_keys:
        if other_key not in mut_inv_map:
            mut_index_map[current_idx] = other_key
            mut_inv_map[other_key] = current_idx
            current_idx += 1

    # Compute distances: each mutation -> all other nodes (other mutations and residue nodes)
    all_mut_edges: List[Tuple[int, int, float]] = []
    for mut_idx, mut_key in enumerate(mut_keys_list):
        mut_atoms = res_atoms_all.get(mut_key)
        if not mut_atoms:
            continue

        # Connect to all other nodes (other mutations and residue nodes)
        all_other_keys = set(mut_keys_list) | set(other_residue_keys)
        all_other_keys.discard(mut_key)  # Exclude self

        for other_key in sorted(all_other_keys):
            other_atoms = res_atoms_all.get(other_key)
            if not other_atoms:
                # Log warning if this is a mutation node (indices 0, 1, 2, ...)
                # This indicates a bug where the mutation key doesn't match what's in res_atoms_all
                if other_key in mut_keys_list:
                    failure_logger.warning(
                        f"Mutation node {other_key} not found in res_atoms_all when computing distance from "
                        f"mutation node {mut_key}. This may indicate a key normalization bug. "
                        f"Available keys in res_atoms_all for this position: "
                        f"{[k for k in res_atoms_all.keys() if str(k[0]).strip().upper() == str(other_key[0]).strip().upper() and normalize_resseq(str(k[1]).strip()) == normalize_resseq(str(other_key[1]).strip())]}"
                    )
                continue

            dval = residue_distance(mut_atoms, other_atoms, contact_distance_strategy)
            other_idx = mut_inv_map[other_key]
            all_mut_edges.append((mut_idx, other_idx, dval))

    # Log completion with summary statistics
    num_mutations = len(mut_keys_list)
    num_edges = len(all_mut_edges)
    logger.info(
        f"Finished mutation edge calculation: {num_mutations} mutation(s) processed, " f"{num_edges} edge(s) computed"
    )

    return mut_index_map, all_mut_edges
