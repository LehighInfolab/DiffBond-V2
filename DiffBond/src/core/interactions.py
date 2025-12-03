"""
interactions.py
===============

Core interaction modes for DiffBond.

This module implements functions to detect specific types of molecular
interactions between two protein chains or structures. Each mode corresponds
to a different class of interaction, using geometric and/or electrostatic
criteria:

    - Contact (c_mode): generic atom contacts within a distance cutoff
    - Ionic (i_mode): electrostatically defined ionic bonds
    - Adjacent (a_mode): atoms adjacent to an existing bond set
    - Cation–π (p_mode): interactions between cations and aromatic rings
    - Salt bridge (s_mode): combined ionic + hydrogen bond pairs
    - Hydrogen bond (h_mode): pre-filter then HBondFinder analysis

All modes return a consistent tuple structure
    (atom_index, edge_list, raw_edges)
which can be used downstream to build graphs or perform further analysis.
"""

import os
from pathlib import Path
import shutil
import tempfile
import logging
from typing import List, Tuple, Optional, Any, Dict

import src.utils.file_handler as file_handler
import src.core.hbondfinder_handler as hbondfinder_handler
import src.utils.distance as distance
from src.utils.graph_handler import convert_to_indexed_edge_list

logger = logging.getLogger(__name__)


def _residue_key(atom: List[str]) -> Tuple[str, str, str]:
    """Return a unique residue identifier (chain, resSeq, resName)."""
    # atom indices: [.., resName(3), atomName(2), chain(4), resSeq(5), x(6), y(7), z(8) ...]
    return atom[4], atom[5], atom[3]


def _build_residue_atoms(points: List[List[str]]) -> Dict[Tuple[str, str, str], List[List[str]]]:
    by_res: Dict[Tuple[str, str, str], List[List[str]]] = {}
    for a in points:
        key = _residue_key(a)
        by_res.setdefault(key, []).append(a)
    return by_res


def _centroid(coords: List[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    import numpy as np

    arr = np.array(coords, dtype=float)
    c = arr.mean(axis=0)
    return float(c[0]), float(c[1]), float(c[2])


def _coords(atom: List[str]) -> Tuple[float, float, float]:
    return float(atom[6]), float(atom[7]), float(atom[8])


def c_mode(
    PDB_data: List[Any],
    dist: float,
    verbose: bool = False,
    node_level: str = "atom",
    contact_distance_strategy: str = "avg-atom",
):
    """
    Detect generic atom contacts within a cutoff distance.

    Args:
        PDB_data (list): Two parsed PDB chain/dataset objects: [entity_A, entity_B].
        dist (float): Maximum distance in Å to consider an atom pair a contact.
        verbose (bool): If True, log detailed edge lists.

    Returns:
        tuple: (atom_index, edge_list, contact_edges)
            atom_index (dict): atom → integer index
            edge_list (list[tuple[int,int]]): edges using integer indices
            contact_edges (list[tuple]): raw atom pair data
    """
    logger.info("Searching contacts within %.2f Å...", dist)
    contact_edges = distance.find_electrostatic_interactions(PDB_data[0], PDB_data[1], dist, False)

    if node_level == "atom":
        atom_index, edge_list = convert_to_indexed_edge_list(contact_edges)
        # augment edge list with distances
        inv_map: Dict[Tuple, int] = {tuple(atom): idx for idx, atom in atom_index.items()}
        edge_list_with_d = []
        for a1, a2 in contact_edges:
            i = inv_map[tuple(a1)]
            j = inv_map[tuple(a2)]
            d = distance.euclidean_distance(a1, a2)
            edge_list_with_d.append((i, j, d))
        if verbose:
            logger.debug("--- CONTACT edges (with distances) ---\n%s", edge_list_with_d)
        return atom_index, edge_list_with_d, contact_edges

    # residue-level aggregation
    res_atoms1 = _build_residue_atoms(PDB_data[0])
    res_atoms2 = _build_residue_atoms(PDB_data[1])

    # group atom pairs by residue pair
    grouped: Dict[Tuple[Tuple[str, str, str], Tuple[str, str, str]], List[Tuple[List[str], List[str]]]] = {}
    for a1, a2 in contact_edges:
        k = (_residue_key(a1), _residue_key(a2))
        grouped.setdefault(k, []).append((a1, a2))

    # compute residue-residue distance per strategy
    residue_edges: List[Tuple[Tuple[str, str, str], Tuple[str, str, str], float]] = []
    for (r1, r2), pairs in grouped.items():
        if contact_distance_strategy == "avg-atom":
            dvals = [distance.euclidean_distance(a1, a2) for (a1, a2) in pairs]
            dval = sum(dvals) / len(dvals)
        else:  # centroid
            c1 = _centroid([_coords(a) for a in res_atoms1[r1]])
            c2 = _centroid([_coords(a) for a in res_atoms2[r2]])
            # Build pseudo-atoms with centroid coordinates in the expected format positions
            # Just reuse first atom as template for euclidean distance function by injecting coords
            # but easier: compute directly
            from math import sqrt

            dval = sqrt((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2)
        residue_edges.append((r1, r2, dval))

    # index residues
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
    edge_list_with_d = [(inv_map[r1], inv_map[r2], d) for (r1, r2, d) in residue_edges]

    if verbose:
        logger.debug("--- CONTACT residue edges (with distances) ---\n%s", edge_list_with_d)
    return res_index, edge_list_with_d, residue_edges


def i_mode(
    PDB_data: List[Any],
    dist: float,
    verbose: bool = False,
    node_level: str = "atom",
):
    """
    Detect ionic interactions using electrostatic criteria.

    Args:
        PDB_data (list): Two parsed PDB chain/dataset objects.
        dist (float): Maximum distance in Å.
        verbose (bool): If True, log atom indices and edge lists.

    Returns:
        tuple: (atom_index, edge_list, ionic_edges)
    """
    logger.info("Searching ionic bonds (%.2f Å cutoff)...", dist)
    ionic_edges = distance.find_electrostatic_interactions(PDB_data[0], PDB_data[1], dist, True)

    if node_level == "atom":
        atom_index, edge_list = convert_to_indexed_edge_list(ionic_edges)
        if verbose:
            logger.debug("--- IONIC atom index ---\n%s", atom_index)
            logger.debug("--- IONIC edges ---\n%s", edge_list)
        return atom_index, edge_list, ionic_edges

    # residue-level: choose distance from actual relevant atom pair (min over pairs)
    grouped: Dict[Tuple[Tuple[str, str, str], Tuple[str, str, str]], List[Tuple[List[str], List[str]]]] = {}
    for a1, a2 in ionic_edges:
        k = (_residue_key(a1), _residue_key(a2))
        grouped.setdefault(k, []).append((a1, a2))

    residue_edges: List[Tuple[Tuple[str, str, str], Tuple[str, str, str], float]] = []
    for (r1, r2), pairs in grouped.items():
        dmin = min(distance.euclidean_distance(a1, a2) for (a1, a2) in pairs)
        residue_edges.append((r1, r2, dmin))

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
    edge_list_with_d = [(inv_map[r1], inv_map[r2], d) for (r1, r2, d) in residue_edges]

    if verbose:
        logger.debug("--- IONIC residue edges (min atom distance) ---\n%s", edge_list_with_d)
    return res_index, edge_list_with_d, residue_edges


def a_mode(PDB_data: List[Any], dist: float, bonds: List[Any], verbose: bool = False):
    """
    Detect atoms adjacent to an existing set of bonds.

    Args:
        PDB_data (list): Two parsed PDB chain/dataset objects.
        dist (float): Cutoff distance in Å.
        bonds (list): List of existing edges/bonds to compare against.
        verbose (bool): If True, log details of adjacent edges.

    Returns:
        tuple: (atom_index, edge_list, adj_edges)
    """
    logger.info("Searching contacts adjacent to bonds (%.2f Å)...", dist)
    # Use nearby-contacts utility to expand around the supplied bond endpoints
    adj_edges = distance.find_nearby_contacts(bonds, PDB_data[0], PDB_data[1], dist)

    # Index and attach distances like contact mode (atom-level)
    atom_index, edge_list = convert_to_indexed_edge_list(adj_edges)
    inv_map: Dict[Tuple, int] = {tuple(atom): idx for idx, atom in atom_index.items()}
    edge_list_with_d = []
    for a1, a2 in adj_edges:
        i = inv_map[tuple(a1)]
        j = inv_map[tuple(a2)]
        d = distance.euclidean_distance(a1, a2)
        edge_list_with_d.append((i, j, d))

    if verbose:
        logger.debug("--- ADJACENT edges (with distances) ---\n%s", edge_list_with_d)
    return atom_index, edge_list_with_d, adj_edges


def p_mode(
    PDB_data: List[Any],
    dist: float,
    ring_radius: float,
    verbose: bool = False,
    use_cylinder: bool = True,
    cylinder_height: float = 6.0,
    node_level: str = "atom",
):
    """
    Detect cation–π interactions between charged groups and aromatic rings.

    Args:
        PDB_data (list): Two parsed PDB chain/dataset objects.
        dist (float): Cutoff distance in Å (or cylinder radius if use_cylinder=True).
        ring_radius (float): Expected aromatic ring radius (or cylinder radius if use_cylinder=True).
        verbose (bool): If True, log details of cation–π interactions.
        use_cylinder (bool): If True, use cylinder-based method (default: True).
        cylinder_height (float): Height of cylinder extending on either side of ring plane (default: 6.0).
        node_level (str): Granularity level ('atom' or 'residue', default: 'atom').

    Returns:
        tuple: (atom_index, edge_list, cation_pi_edges)
    """
    logger.info("Searching cation–pi interactions...")

    if use_cylinder:
        # Use cylinder-based method
        # ring_radius parameter is used as cylinder_radius
        cation_pi_edges = distance.find_cation_pi_interactions_cylinder(
            PDB_data[0],
            PDB_data[1],
            cylinder_radius=ring_radius,
            cylinder_height=cylinder_height,
        )
    else:
        # Use original distance-based method
        cation_pi_edges = distance.find_cation_pi_interactions(PDB_data[0], PDB_data[1], dist + 1, ring_radius)

    if node_level == "atom":
        atom_index, edge_list = convert_to_indexed_edge_list(cation_pi_edges)

        # Augment edge list with distances
        inv_map: Dict[Tuple, int] = {tuple(atom): idx for idx, atom in atom_index.items()}
        edge_list_with_d = []
        for a1, a2 in cation_pi_edges:
            i = inv_map.get(tuple(a1))
            j = inv_map.get(tuple(a2))
            if i is not None and j is not None:
                d = distance.euclidean_distance(a1, a2)
                edge_list_with_d.append((i, j, d))

        if verbose:
            logger.debug("--- CATION-π atom index ---\n%s", atom_index)
            logger.debug("--- CATION-π edges (with distances) ---\n%s", edge_list_with_d)

        return atom_index, edge_list_with_d, cation_pi_edges

    # residue-level: choose distance from actual relevant atom pair (min over pairs)
    grouped: Dict[Tuple[Tuple[str, str, str], Tuple[str, str, str]], List[Tuple[List[str], List[str]]]] = {}
    for a1, a2 in cation_pi_edges:
        k = (_residue_key(a1), _residue_key(a2))
        grouped.setdefault(k, []).append((a1, a2))

    residue_edges: List[Tuple[Tuple[str, str, str], Tuple[str, str, str], float]] = []
    for (r1, r2), pairs in grouped.items():
        dmin = min(distance.euclidean_distance(a1, a2) for (a1, a2) in pairs)
        residue_edges.append((r1, r2, dmin))

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
    edge_list_with_d = [(inv_map[r1], inv_map[r2], d) for (r1, r2, d) in residue_edges]

    if verbose:
        logger.debug("--- CATION-π residue edges (min atom distance) ---\n%s", edge_list_with_d)
    return res_index, edge_list_with_d, residue_edges


def s_mode(ionic_edges: List[Any], hbond_edges: List[Any], verbose: bool = False):
    """
    Detect salt bridges as overlapping ionic and hydrogen-bond interactions.

    Args:
        ionic_edges (list): Previously detected ionic interactions.
        hbond_edges (list): Previously detected hydrogen bonds.
        verbose (bool): If True, log details of matched salt bridge edges.

    Returns:
        tuple: ([atom_index1, atom_index2],
                [edge_list1, edge_list2],
                [saltbridge_i_edges, saltbridge_h_edges])
    """
    if ionic_edges is None or hbond_edges is None:
        logger.warning("Salt bridge detection requires both ionic and hbond edges.")
        return None, None, None

    def is_residue_edge_format(item: Any) -> bool:
        # Residue edges we produce are (resKey, resKey, distance)
        return isinstance(item, (list, tuple)) and len(item) == 3 and isinstance(item[0], (list, tuple))

    if (
        ionic_edges
        and is_residue_edge_format(ionic_edges[0])
        and hbond_edges
        and is_residue_edge_format(hbond_edges[0])
    ):
        # Residue-level: intersect residue pairs
        def pair_key(e):
            return (e[0], e[1])

        ionic_map = {pair_key(e): e for e in ionic_edges}
        hbond_map = {pair_key(e): e for e in hbond_edges}

        common_pairs = set(ionic_map.keys()) & set(hbond_map.keys())
        saltbridge_i_edges = [ionic_map[k] for k in common_pairs]
        saltbridge_h_edges = [hbond_map[k] for k in common_pairs]

        # Build residue indices
        res_index: Dict[int, Tuple[str, str, str]] = {}
        current = 0
        for r1, r2, _ in saltbridge_i_edges + saltbridge_h_edges:
            if r1 not in res_index.values():
                res_index[current] = r1
                current += 1
            if r2 not in res_index.values():
                res_index[current] = r2
                current += 1

        inv_map: Dict[Tuple[str, str, str], int] = {v: k for k, v in res_index.items()}
        edge_list1 = [(inv_map[r1], inv_map[r2], d) for (r1, r2, d) in saltbridge_i_edges]
        edge_list2 = [(inv_map[r1], inv_map[r2], d) for (r1, r2, d) in saltbridge_h_edges]

        if verbose:
            logger.debug("--- SALT BRIDGE (residue) ionic edges ---\n%s", edge_list1)
            logger.debug("--- SALT BRIDGE (residue) hbond edges ---\n%s", edge_list2)

        return [res_index, res_index], [edge_list1, edge_list2], [saltbridge_i_edges, saltbridge_h_edges]

    # Atom-level fallback: use previous matching by residue number fields at [5]
    def match_amino_acid_tuples(list1, list2):
        """Match residue pairs that overlap in ionic and hbond lists (atom-level inputs)."""
        l1_matches, l2_matches = [], []
        for pair1 in list1:
            for pair2 in list2:
                paired = (
                    pair1[0][5].strip() == pair2[0][5].strip() and pair1[1][5].strip() == pair2[1][5].strip()
                ) or (pair1[0][5].strip() == pair2[1][5].strip() and pair1[1][5].strip() == pair2[0][5].strip())
                if paired:
                    if pair1 not in l1_matches:
                        l1_matches.append(pair1)
                    if pair2 not in l2_matches:
                        l2_matches.append(pair2)
        return l1_matches, l2_matches

    saltbridge_i_edges, saltbridge_h_edges = match_amino_acid_tuples(ionic_edges, hbond_edges)

    atom_index1, edge_list1 = convert_to_indexed_edge_list(saltbridge_i_edges)
    atom_index2, edge_list2 = convert_to_indexed_edge_list(saltbridge_h_edges)

    if verbose:
        logger.debug("--- SALT BRIDGE ionic edges ---\n%s", edge_list1)
        logger.debug("--- SALT BRIDGE hbond edges ---\n%s", edge_list2)

    return [atom_index1, atom_index2], [edge_list1, edge_list2], [saltbridge_i_edges, saltbridge_h_edges]


def cov_mode(
    PDB_data: List[Any],
    allowed_residue_keys: set,
    verbose: bool = False,
    node_level: str = "residue",
):
    """Build covalent backbone edges between consecutive residues (i,i+1)
    restricted to residues that appear in the contact node set.

    Args:
        PDB_data: Two parsed PDB chain/dataset objects (list of atoms per entity).
        allowed_residue_keys: Set of residue keys (chain, resSeq, resName) to include as nodes.
        verbose: Log details if True.
        node_level: Only "residue" is supported; atom-level is not applicable here.

    Returns:
        tuple: (res_index, edge_list_with_d, residue_edges)
            - res_index: dict[int -> (chain,resSeq,resName)] over allowed residues encountered
            - edge_list_with_d: list[(i,j,d)] edges between consecutive residues within a chain
              where d is computed using contact_distance_strategy
            - residue_edges: list[((chain,resSeq,resName),(chain,resSeq,resName),d)]
    """

    # Normalize keys for robust matching
    def _norm_key(rk):
        return (str(rk[0]).strip(), str(rk[1]).strip(), str(rk[2]).strip())

    allowed = {_norm_key(k) for k in allowed_residue_keys}

    # Collect residues per chain from both entities, restricted to allowed set
    by_chain: Dict[str, set] = {}

    def _collect(points: List[List[str]]):
        for a in points:
            key = _norm_key((a[4], a[5], a[3]))
            if key in allowed:
                by_chain.setdefault(key[0], set()).add(key)

    _collect(PDB_data[0])
    _collect(PDB_data[1])

    # Build residue->atoms map across both entities for distance computation
    res_atoms_all: Dict[Tuple[str, str, str], List[List[str]]] = {}

    def _add_residue_atoms(points: List[List[str]]):
        for a in points:
            key = _norm_key((a[4], a[5], a[3]))
            res_atoms_all.setdefault(key, []).append(a)

    _add_residue_atoms(PDB_data[0])
    _add_residue_atoms(PDB_data[1])

    # Distance between two consecutive residues: use C(i)–N(i+1) peptide bond distance
    from math import sqrt as _sqrt

    def _centroid_of(atoms: List[List[str]]) -> Tuple[float, float, float]:
        xs = [float(a[6]) for a in atoms]
        ys = [float(a[7]) for a in atoms]
        zs = [float(a[8]) for a in atoms]
        return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)

    def _atom_coord(atoms: List[List[str]], atom_name: str) -> Optional[Tuple[float, float, float]]:
        for a in atoms:
            try:
                if str(a[2]).strip().upper() == atom_name:
                    return float(a[6]), float(a[7]), float(a[8])
            except Exception:
                continue
        return None

    def _residue_distance(a_list: List[List[str]], b_list: List[List[str]]) -> float:
        c_coord = _atom_coord(a_list, "C")
        n_coord = _atom_coord(b_list, "N")
        if c_coord is not None and n_coord is not None:
            dx = c_coord[0] - n_coord[0]
            dy = c_coord[1] - n_coord[1]
            dz = c_coord[2] - n_coord[2]
            return _sqrt(dx * dx + dy * dy + dz * dz)
        # Fallbacks if specific atoms are missing
        ca1 = _atom_coord(a_list, "CA")
        ca2 = _atom_coord(b_list, "CA")
        if ca1 is not None and ca2 is not None:
            return _sqrt((ca1[0] - ca2[0]) ** 2 + (ca1[1] - ca2[1]) ** 2 + (ca1[2] - ca2[2]) ** 2)
        c1 = _centroid_of(a_list)
        c2 = _centroid_of(b_list)
        return _sqrt((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2)

    # Sort residues within each chain by residue sequence (as integer if possible)
    def _seq_num(rk):
        try:
            return int(str(rk[1]).strip())
        except Exception:
            return str(rk[1]).strip()

    residue_edges: List[Tuple[Tuple[str, str, str], Tuple[str, str, str], float]] = []
    for chain, res_set in by_chain.items():
        ordered = sorted(res_set, key=_seq_num)
        for i in range(len(ordered) - 1):
            r1 = ordered[i]
            r2 = ordered[i + 1]
            # Ensure residue numbers are actually consecutive (i and i+1)
            try:
                seq1 = int(str(r1[1]).strip())
                seq2 = int(str(r2[1]).strip())
                if seq2 != seq1 + 1:
                    continue  # Skip if not consecutive
            except (ValueError, IndexError):
                # If we can't parse as integers, skip this pair
                continue
            a_list = res_atoms_all.get(r1)
            b_list = res_atoms_all.get(r2)
            if not a_list or not b_list:
                continue
            dval = _residue_distance(a_list, b_list)
            residue_edges.append((r1, r2, dval))

    # Build residue index and convert to (i,j,d)
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

    edge_list_with_d = [(inv_map[r1], inv_map[r2], d) for (r1, r2, d) in residue_edges]

    if verbose:
        logger.debug("--- COVALENT residue edges (residue distances) ---\n%s", edge_list_with_d)
    return res_index, edge_list_with_d, residue_edges


def h_mode(
    PDB_data: List[Any],
    outputFileName: str,
    verbose: bool = False,
    node_level: str = "atom",
):
    """
    Detect hydrogen bonds using a prefilter and HBondFinder.

    Args:
        PDB_data (list): Two parsed PDB chain/dataset objects.
        outputFileName (str): Name for the output subdirectory under Results/.
        verbose (bool): If True, log details of hydrogen bonds.

    Returns:
        tuple: (atom_index, edge_list, hbond_edges)
    """
    logger.info("Searching hydrogen bonds (prefilter 3.5 Å)...")
    edges_temp = distance.find_electrostatic_interactions(PDB_data[0], PDB_data[1], 3.5, False)

    logger.info("Running HBondFinder...")
    with tempfile.TemporaryDirectory(delete=False) as temp_dir:
        HB_file, hb_file = handle_hbond_processing(edges_temp, Path(temp_dir), outputFileName)

    if HB_file is None or hb_file is None:
        logger.error("No eligible contacts to run HBondFinder.")
        return {}, [], []

    atom_index, edge_list = hbondfinder_handler.parse_hblines_file(hb_file)
    hbond_edges = hbondfinder_handler.edges_with_atom_details(atom_index, edge_list)

    if node_level == "atom":
        if verbose:
            logger.debug("--- H-BOND edges ---\n%s", hbond_edges)
        return atom_index, edge_list, hbond_edges

    # residue-level: aggregate bonds by residue pair, distance from actual donor/acceptor atoms
    from math import sqrt

    def euclid_from_atoms(a1: List[str], a2: List[str]) -> float:
        x1, y1, z1 = float(a1[6]), float(a1[7]), float(a1[8])
        x2, y2, z2 = float(a2[6]), float(a2[7]), float(a2[8])
        return sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)

    grouped: Dict[Tuple[Tuple[str, str, str], Tuple[str, str, str]], List[Tuple[List[str], List[str]]]] = {}
    for a1, a2 in hbond_edges:
        k = (_residue_key(a1), _residue_key(a2))
        grouped.setdefault(k, []).append((a1, a2))

    residue_edges: List[Tuple[Tuple[str, str, str], Tuple[str, str, str], float]] = []
    for (r1, r2), pairs in grouped.items():
        dmin = min(euclid_from_atoms(a1, a2) for (a1, a2) in pairs)
        residue_edges.append((r1, r2, dmin))

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
    edge_list_with_d = [(inv_map[r1], inv_map[r2], d) for (r1, r2, d) in residue_edges]

    if verbose:
        logger.debug("--- H-BOND residue edges (min atom distance) ---\n%s", edge_list_with_d)
    return res_index, edge_list_with_d, residue_edges


def handle_hbond_processing(
    edges: List[Any], temp_dir: Path, output_file_name: str
) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Prepare temporary PDBs, run HBondFinder, and move outputs to Results/.

    Args:
        edges (list): Candidate contact edges from prefilter step.
        temp_dir (Path): Path to a temporary working directory.
        output_file_name (str): Subdirectory name under Results/.

    Returns:
        tuple: (HBondFinder_output.txt path, hbonds_list.txt path), or (None, None)
    """
    if not edges:
        return None, None

    atoms1 = [edge[0] for edge in edges]
    atoms2 = [edge[1] for edge in edges]

    # Write atoms to PDB files
    temp_pdb_path = temp_dir / "temp.pdb"
    file_handler.write_PDB(str(temp_pdb_path), False, atoms1)
    file_handler.write_PDB(str(temp_pdb_path), True, atoms2)

    output_dir = Path("results") / output_file_name

    # Run hydrogen bond finder
    try:
        hbondfinder_handler.run_hbondfinder(str(temp_pdb_path))
    except Exception as error:
        logger.exception(f"Failed to run hbondfinder {error}")
        output_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(temp_pdb_path, output_dir / "hbondfinder_error.txt")
        return None, None

    try:
        # Move HBondFinder files
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.debug("Attempting to move HBondFinder_temp.txt...")
        shutil.copy("HBondFinder_temp.txt", output_dir / "HBondFinder_output.txt")
        shutil.copy("hbonds_temp.txt", output_dir / "hbonds_list.txt")
    except OSError as error:
        logger.exception(f"Failed to copy HBondFinder output files {error}")
        return None, None

    try:
        logger.debug("Attempting to remove temp files...")
        Path.unlink("HBondFinder_temp.txt")
        Path.unlink("hbonds_temp.txt")
    except Exception as error:
        logger.warning(f"Failed to remove HBondFinder output files {error}")

    return output_dir / "HbondFinder_output.txt", output_dir / "hbonds_list.txt"


def handle_interchain_hbond_processing(file: str) -> Optional[str]:
    """
    Legacy: run HBondFinder on a single PDB file and move outputs.

    Args:
        file (str): Path to a PDB file.

    Returns:
        str: Name of the HBondFinder output file, or None if failed.
    """
    try:
        shutil.copyfile(file, "./temp.pdb")
    except OSError:
        logger.error("File not found: %s", file)
        return None

    if not hbondfinder_handler.run_hbondfinder("temp.pdb"):
        os.remove("temp.pdb")
        return None

    hb_file_name = f"{file.split('.')[-2].split('\\')[-1]}.txt"
    try:
        os.rename("HBondFinder_temp.txt", f"HBondFinder{hb_file_name}")
        shutil.move(f"HBondFinder{hb_file_name}", "hbondfinder_data")
        os.rename("hbonds_temp.txt", f"hbonds{hb_file_name}")
        shutil.move(f"hbonds{hb_file_name}", "hbondfinder_data")
    except OSError:
        logger.error("Could not move HBondFinder outputs to hbondfinder_data.")
    finally:
        os.remove("temp.pdb")

    return f"HBondFinder{hb_file_name}"


def compute_intra_contact_edges(
    pdb_data: List[Any],
    contact_index: Dict[int, List[Any]] | Dict[int, Tuple[str, str, str]],
    contact_raw: List[Any],
    max_distance: float,
    node_level: str,
    contact_distance_strategy: str,
) -> List[Tuple[int, int, float]]:
    """Compute intramolecular contact edges restricted to nodes present in contact edges.

    Returns:
        List of edge triplets (source_index, target_index, distance).
    """
    if node_level == "atom":
        atoms_in_contacts = set()
        inv_map: Dict[Tuple, int] = {tuple(atom): idx for idx, atom in contact_index.items()}

        for a1, a2 in contact_raw:
            idx1 = inv_map.get(tuple(a1))
            idx2 = inv_map.get(tuple(a2))
            if idx1 is not None:
                atoms_in_contacts.add(idx1)
            if idx2 is not None:
                atoms_in_contacts.add(idx2)

        def _get_chain_from_idx(idx: int) -> str:
            atom = contact_index[idx]
            return atom[4].strip() if len(atom) > 4 else ""

        chain1 = _get_chain_from_idx(0) if 0 in contact_index else ""
        chain2 = None
        if len(contact_index) > 0:
            for idx in contact_index:
                ch = _get_chain_from_idx(idx)
                if ch != chain1:
                    chain2 = ch
                    break

        ent1_indices: List[int] = []
        ent2_indices: List[int] = []
        for idx in atoms_in_contacts:
            ch = _get_chain_from_idx(idx)
            if ch == chain1:
                ent1_indices.append(idx)
            elif chain2 is not None and ch == chain2:
                ent2_indices.append(idx)

        def _intra_pairs(indices: List[int]) -> List[Tuple[int, int, float]]:
            res: List[Tuple[int, int, float]] = []
            atoms_list = [contact_index[i] for i in indices]
            for i in range(len(atoms_list)):
                for j in range(i + 1, len(atoms_list)):
                    d = distance.euclidean_distance(atoms_list[i], atoms_list[j])
                    if d < max_distance:
                        res.append((indices[i], indices[j], d))
            return res

        intra_edges = _intra_pairs(ent1_indices) + _intra_pairs(ent2_indices)
        return intra_edges

    # residue-level
    def _normalize_res_key(rk: Tuple[str, str, str]) -> Tuple[str, str, str]:
        return (str(rk[0]).strip(), str(rk[1]).strip(), str(rk[2]).strip())

    union_res: set = set()
    for r1, r2, _ in contact_raw:
        union_res.add(_normalize_res_key(r1))
        union_res.add(_normalize_res_key(r2))

    def _build_residue_atoms(points: List[List[str]]):
        by_res: Dict[Tuple[str, str, str], List[List[str]]] = {}
        for a in points:
            key = _normalize_res_key((a[4], a[5], a[3]))
            by_res.setdefault(key, []).append(a)
        return by_res

    res_atoms: Dict[Tuple[str, str, str], List[List[str]]] = {}
    res_atoms.update(_build_residue_atoms(pdb_data[0]))
    res_atoms.update(_build_residue_atoms(pdb_data[1]))

    inv_idx: Dict[Tuple[str, str, str], int] = {}
    for k, v in contact_index.items():
        norm_key = _normalize_res_key(v)
        inv_idx[norm_key] = k

    from math import sqrt as _sqrt2

    def _centroid_of2(atoms: List[List[str]]):
        xs = [float(a[6]) for a in atoms]
        ys = [float(a[7]) for a in atoms]
        zs = [float(a[8]) for a in atoms]
        return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)

    def _residue_distance(a_list: List[List[str]], b_list: List[List[str]]):
        if contact_distance_strategy == "avg-atom":
            vals: List[float] = []
            for a1 in a_list:
                x1, y1, z1 = float(a1[6]), float(a1[7]), float(a1[8])
                for a2 in b_list:
                    x2, y2, z2 = float(a2[6]), float(a2[7]), float(a2[8])
                    vals.append(_sqrt2((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2))
            return (sum(vals) / len(vals)) if vals else 0.0
        c1 = _centroid_of2(a_list)
        c2 = _centroid_of2(b_list)
        return _sqrt2((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2)

    by_chain: Dict[str, List[Tuple[str, str, str]]] = {}
    for rk in union_res:
        by_chain.setdefault(rk[0], []).append(rk)

    intra_edges: List[Tuple[int, int, float]] = []
    for chain, keys in by_chain.items():
        keys_sorted = sorted(keys)
        for i in range(len(keys_sorted)):
            for j in range(i + 1, len(keys_sorted)):
                a_list = res_atoms.get(keys_sorted[i])
                b_list = res_atoms.get(keys_sorted[j])
                if not a_list or not b_list:
                    continue
                d = _residue_distance(a_list, b_list)
                if d < max_distance:
                    ii = inv_idx.get(keys_sorted[i])
                    jj = inv_idx.get(keys_sorted[j])
                    if ii is not None and jj is not None:
                        intra_edges.append((ii, jj, d))

    return intra_edges
