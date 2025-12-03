"""
Utility functions for processing SKEMPI mutations and computing mutant-to-residue distances.
"""

import csv
import re
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Set, Any
from math import sqrt

from src.core.interactions import _residue_key, _build_residue_atoms

logger = logging.getLogger(__name__)


def parse_mutation_token(token: str) -> Tuple[str, str] | None:
    """Parse a mutation token like 'LI45G' into (chain, resSeq).

    Args:
        token: Mutation token in format like 'LI45G'

    Returns:
        Tuple of (chain, resSeq) or None if parsing fails
    """
    m = re.match(r"^[A-Z]([A-Z])(\d+)[A-Z]$", token)
    if not m:
        return None
    return m.group(1), m.group(2)


def read_skempi_mutations(skempi_csv: Path, row_index: int) -> List[str]:
    """Read mutation specifications from a SKEMPI CSV file.

    Args:
        skempi_csv: Path to SKEMPI CSV file (semicolon-delimited)
        row_index: Row index to read (header skipped)

    Returns:
        List of mutation token strings (may be comma-separated in CSV, will be split)
    """
    mut_specs: List[str] = []
    try:
        with open(skempi_csv, "r") as f:
            reader = csv.reader(f, delimiter=";")
            for i, row in enumerate(reader):
                if i == 0:
                    continue  # header
                if i == row_index:
                    # Typical SKEMPI v2 columns: 0=PDB id, 1=mutant on chain1, 2=mutant on chain2
                    # if len(row) > 1 and row[1].strip():
                    #     mut_specs.append(row[1].strip())
                    if len(row) > 2 and row[2].strip():
                        # Split by comma to handle multiple mutations
                        tokens = [
                            t.strip() for t in row[2].strip().split(",") if t.strip()
                        ]
                        mut_specs.extend(tokens)
                    break
    except Exception as err:
        logger.exception("Failed reading SKEMPI CSV for mutant mapping: %s", err)
        mut_specs = []
    return mut_specs


def centroid_of(atoms: List[List[str]]) -> Tuple[float, float, float]:
    """Calculate centroid of a list of atoms.

    Args:
        atoms: List of atom data (PDB format)

    Returns:
        Tuple of (x, y, z) coordinates
    """
    xs = [float(a[6]) for a in atoms]
    ys = [float(a[7]) for a in atoms]
    zs = [float(a[8]) for a in atoms]
    return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)


def residue_distance(
    mut_atoms: List[List[str]],
    other_atoms: List[List[str]],
    strategy: str = "avg-atom",
) -> float:
    """Calculate distance between two residues.

    Args:
        mut_atoms: Atoms of first residue
        other_atoms: Atoms of second residue
        strategy: 'avg-atom' for average atom-atom distance, 'centroid' for centroid distance

    Returns:
        Distance in Angstroms
    """
    if strategy == "avg-atom":
        dvals = []
        for a1 in mut_atoms:
            x1, y1, z1 = float(a1[6]), float(a1[7]), float(a1[8])
            for a2 in other_atoms:
                x2, y2, z2 = float(a2[6]), float(a2[7]), float(a2[8])
                dvals.append(sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2))
        return (sum(dvals) / len(dvals)) if dvals else 0.0

    # centroid strategy
    c1 = centroid_of(mut_atoms)
    c2 = centroid_of(other_atoms)
    return sqrt((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2)


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
                    keys.add(_residue_key(a1))
                if isinstance(a2, list):
                    keys.add(_residue_key(a2))

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
        union_residue_keys.update(
            _extract_residue_keys_from_edges(edge_dict["ionic"][2])
        )

    # From hbond (if present)
    if "hbond" in edge_dict and edge_dict["hbond"][2]:
        union_residue_keys.update(
            _extract_residue_keys_from_edges(edge_dict["hbond"][2])
        )

    # From salt bridges
    if "saltbridge_ionic" in edge_dict or "saltbridge_hbond" in edge_dict:
        for sb_type in ["saltbridge_ionic", "saltbridge_hbond"]:
            if sb_type in edge_dict:
                sb_raw = edge_dict[sb_type][2]
                if sb_raw:
                    union_residue_keys.update(_extract_residue_keys_from_edges(sb_raw))

    # From adjacent (future mode), if present
    if "adjacent" in edge_dict and edge_dict["adjacent"][2]:
        union_residue_keys.update(
            _extract_residue_keys_from_edges(edge_dict["adjacent"][2])
        )

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
    # Build residue-to-atoms mapping for all entities
    res_atoms_all: Dict[Tuple[str, str, str], List[List[str]]] = {}
    for entity in pdb_data:
        res_atoms_all.update(_build_residue_atoms(entity))

    # Build reverse lookup by (chain, resSeq) ignoring resName
    chain_res_to_key: Dict[Tuple[str, str], Tuple[str, str, str]] = {}
    for key in res_atoms_all.keys():
        chain_res_to_key[(str(key[0]).strip(), str(key[1]).strip())] = key

    # First pass: parse and validate all mutations, collect their keys
    mut_keys_list: List[Tuple[str, str, str]] = []
    for token in mut_specs:
        parsed = parse_mutation_token(token)
        if not parsed:
            logger.warning(
                f"Failed to parse mutation token '{token}'. "
                f"Expected format: single letter + chain + residue number + single letter (e.g., 'LI45G'). "
                f"Skipping this mutation."
            )
            continue

        chain, resseq = parsed
        mut_key = chain_res_to_key.get((chain, resseq))
        if mut_key is None:
            logger.warning(
                f"Could not find mutant residue for token '{token}' (chain={chain}, resSeq={resseq}) "
                f"in PDB structure. The residue may not exist in the input PDB files."
            )
            continue

        mut_atoms = res_atoms_all.get(mut_key)
        if not mut_atoms:
            logger.warning(
                f"No atoms found for mutant residue {mut_key} from token '{token}'. "
                f"This may indicate a data parsing issue."
            )
            continue

        mut_keys_list.append(mut_key)

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
                continue

            dval = residue_distance(mut_atoms, other_atoms, contact_distance_strategy)
            other_idx = mut_inv_map[other_key]
            all_mut_edges.append((mut_idx, other_idx, dval))

    return mut_index_map, all_mut_edges
