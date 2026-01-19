"""
Utility functions for handling residue indexing and residue operations.

Main functions:
- normalize_resseq: Normalize residue sequence numbers by removing alternate location suffixes
- position_key: Extract position key (chain, resSeq) from residue key
- build_residue_index: Build residue index including ALL variants as separate nodes
- residue_key: Extract residue identifier from atom data
- build_residue_atoms: Group atoms by residue
"""

import re
from typing import Dict, List, Tuple


def normalize_resseq(resseq: str) -> str:
    """Normalize residue sequence number by removing alternate location suffixes.

    PDB files can have alternate locations like 100, 100A, 100B, 100C.
    This function extracts just the numeric part (e.g., '100A' -> '100').

    Args:
        resseq: Residue sequence number string (may contain suffix like '100A')

    Returns:
        Normalized residue sequence number (numeric part only)
    """
    # Remove any trailing letters (A, B, C, etc.) that indicate alternate locations
    match = re.match(r"^(\d+)", str(resseq).strip())
    if match:
        return match.group(1)
    # If no numeric prefix found, return as-is (shouldn't happen in valid PDB)
    return str(resseq).strip()


def position_key(residue_key: Tuple[str, str, str]) -> Tuple[str, str]:
    """Extract position key (chain, resSeq) from residue key.

    Args:
        residue_key: Tuple of (chain, resSeq, resName)

    Returns:
        Tuple of (chain, resSeq) - position with alternate location suffix preserved
    """
    chain, resseq, _ = residue_key
    return (chain, str(resseq).strip())


def residue_key(atom: List[str]) -> Tuple[str, str, str]:
    """Return a unique residue identifier (chain, resSeq, resName) from atom data.

    Args:
        atom: Atom data list where indices are: [.., resName(3), atomName(2), chain(4), resSeq(5), x(6), y(7), z(8) ...]

    Returns:
        Tuple of (chain, resSeq, resName)
    """
    return atom[4], atom[5], atom[3]


def build_residue_atoms(points: List[List[str]]) -> Dict[Tuple[str, str, str], List[List[str]]]:
    """Group atoms by residue key.

    Args:
        points: List of atom data

    Returns:
        Dictionary mapping residue keys to lists of atoms
    """
    by_res: Dict[Tuple[str, str, str], List[List[str]]] = {}
    for a in points:
        key = residue_key(a)
        by_res.setdefault(key, []).append(a)
    return by_res


def build_residue_index(
    residue_edges: List[Tuple[Tuple[str, str, str], Tuple[str, str, str], float]],
) -> Tuple[
    Dict[int, Tuple[str, str, str]],
    List[Tuple[int, int, float]],
    List[Tuple[Tuple[str, str, str], Tuple[str, str, str], float]],
]:
    """Build residue index including ALL variants (100, 100a, 100b) as separate nodes.

    PDB files can have multiple residues at the same position with insertion codes like 100, 100A, 100B, 100C.
    This function includes ALL variants as separate nodes in the index - no deduplication.
    Each unique (chain, resSeq, resName) combination gets its own index.

    Args:
        residue_edges: List of residue edge tuples (r1, r2, distance)

    Returns:
        Tuple of (res_index, edge_list_with_d, residue_edges)
            - res_index: Dictionary mapping indices to residue keys (all variants included)
            - edge_list_with_d: List of indexed edges with distances
            - residue_edges: List of residue edges (unchanged, all variants preserved)
    """
    # Build index using all unique residue keys (no deduplication)
    res_index: Dict[int, Tuple[str, str, str]] = {}
    inv_map: Dict[Tuple[str, str, str], int] = {}
    current = 0

    for r1, r2, _ in residue_edges:
        if r1 not in inv_map:
            inv_map[r1] = current
            res_index[current] = r1
            current += 1
        if r2 not in inv_map:
            inv_map[r2] = current
            res_index[current] = r2
            current += 1

    # Build edge list with indices
    edge_list_with_d = [(inv_map[r1], inv_map[r2], d) for (r1, r2, d) in residue_edges]

    return res_index, edge_list_with_d, residue_edges
