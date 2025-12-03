"""
Utility functions for calculating distances and finding molecular interactions.

Main functions:
- find_electrostatic_interactions: Find ionic bonds between charged residues
- find_cation_pi_interactions: Find interactions between aromatic rings and cations (distance-based)
- find_cation_pi_interactions_cylinder: Find cation-π interactions using cylinder-based approach
- find_nearby_contacts: Find atoms within specified distance of contact points
"""

from typing import List, Tuple, Dict
import math
import numpy as np
from tqdm import tqdm

from src.utils.constants import (
    AMINO_ACID_PROPERTIES,
    CHARGED_ATOMS,
    AROMATIC_RING_ATOMS,
    CATIONIC_ATOMS,
    PDB_COORDINATE_INDICES as PDB_IDX,
)


def euclidean_distance(point1: List[str], point2: List[str]) -> float:
    """Calculate 3D Euclidean distance between two points.

    Expects points in PDB format where coordinates are at indices defined in PDB_COORDINATE_INDICES.
    """
    return math.sqrt(
        sum((float(point1[PDB_IDX[coord]]) - float(point2[PDB_IDX[coord]])) ** 2 for coord in ["x", "y", "z"])
    )


def find_nearby_contacts(
    contact_points: List[Tuple[List[str], List[str]]],
    points1: List[List[str]],
    points2: List[List[str]],
    max_distance: float,
) -> List[Tuple[List[str], List[str]]]:
    """Find atoms near contact points within specified distance."""
    output = []
    chain1 = points1[0][4]

    for contact in tqdm(contact_points, leave=False):
        # Process contacts where first atom is from chain1
        if contact[0][4] == chain1:
            # Find nearby atoms for first contact atom
            for point in points1:
                if point != contact[0] and euclidean_distance(contact[0], point) < max_distance:
                    output.append((contact[0], point))

            # Find nearby atoms for second contact atom
            for point in points2:
                if point != contact[1] and euclidean_distance(contact[1], point) < max_distance:
                    output.append((contact[1], point))

        # Process contacts where second atom is from chain1
        elif contact[1][4] == chain1:
            for point in points1:
                if point != contact[1] and euclidean_distance(contact[1], point) < max_distance:
                    output.append((contact[1], point))

            for point in points2:
                if point != contact[0] and euclidean_distance(contact[0], point) < max_distance:
                    output.append((contact[0], point))

    return output


def find_electrostatic_interactions(
    points1: List[List[str]],
    points2: List[List[str]],
    max_distance: float,
    check_charge: bool = False,
) -> List[Tuple[List[str], List[str]]]:
    """Find ionic bonds between charged residues."""

    def is_valid_ionic_pair(atom1: List[str], atom2: List[str]) -> bool:
        """Check if two atoms can form an ionic bond."""
        if not check_charge:
            return True

        # Get residue properties
        res1 = AMINO_ACID_PROPERTIES.get(atom1[3], {"charge": 0})
        res2 = AMINO_ACID_PROPERTIES.get(atom2[3], {"charge": 0})

        # Check for opposite charges
        if res1["charge"] * res2["charge"] >= 0:
            return False

        # Verify proper atoms are involved
        atom1_name = atom1[2]
        atom2_name = atom2[2]

        is_pos1 = atom1_name in CHARGED_ATOMS["positive"]
        is_neg1 = atom1_name in CHARGED_ATOMS["negative"]
        is_pos2 = atom2_name in CHARGED_ATOMS["positive"]
        is_neg2 = atom2_name in CHARGED_ATOMS["negative"]

        return (is_pos1 and is_neg2) or (is_neg1 and is_pos2)

    results = []
    for p1 in tqdm(points1, leave=False):
        for p2 in points2:
            if euclidean_distance(p1, p2) < max_distance and is_valid_ionic_pair(p1, p2):
                results.append((p1, p2))

    return results


def find_cation_pi_interactions(
    points1: List[List[str]],
    points2: List[List[str]],
    max_distance: float,
    ring_radius: float,
) -> List[Tuple[List[str], List[str]]]:
    """Find cation-π interactions between aromatic rings and cationic residues."""

    def get_aromatic_center(residue_atoms: List[List[str]]) -> np.ndarray:
        """Calculate center of aromatic ring."""
        coords = [
            [float(atom[PDB_IDX[coord]]) for coord in ["x", "y", "z"]]
            for atom in residue_atoms
            if AMINO_ACID_PROPERTIES.get(atom[3], {}).get("type") == "aromatic"
        ]
        return np.mean(coords, axis=0) if coords else None

    def is_within_volume(point: List[str], center: np.ndarray) -> bool:
        """Check if point is within ring volume."""
        if center is None:
            return False
        coords = [float(point[PDB_IDX[coord]]) for coord in ["x", "y", "z"]]
        return np.linalg.norm(np.array(coords) - center) <= ring_radius

    results = []
    for p1 in tqdm(points1, leave=False):
        if AMINO_ACID_PROPERTIES.get(p1[3], {}).get("type") != "aromatic":
            continue

        aromatic_center = get_aromatic_center(points1)

        for p2 in points2:
            if AMINO_ACID_PROPERTIES.get(p2[3], {}).get("type") != "cationic":
                continue

            if euclidean_distance(p1, p2) < max_distance and is_within_volume(p2, aromatic_center):
                results.append((p1, p2))

    return results


def _residue_key(atom: List[str]) -> Tuple[str, str, str]:
    """Return a unique residue identifier (chain, resSeq, resName)."""
    return atom[4], atom[5], atom[3]


def _get_atom_coords(atom: List[str]) -> np.ndarray:
    """Extract coordinates from atom data."""
    return np.array(
        [
            float(atom[PDB_IDX["x"]]),
            float(atom[PDB_IDX["y"]]),
            float(atom[PDB_IDX["z"]]),
        ]
    )


def _compute_ring_plane(ring_atoms: List[List[str]]) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the centroid and normal vector of an aromatic ring plane.

    Args:
        ring_atoms: List of atoms forming the aromatic ring.

    Returns:
        Tuple of (centroid, normal_vector) where normal_vector is normalized.
    """
    coords = np.array([_get_atom_coords(atom) for atom in ring_atoms])
    centroid = np.mean(coords, axis=0)

    # Compute normal vector using cross product of two ring edges
    # Use first three atoms to define plane
    if len(coords) >= 3:
        v1 = coords[1] - coords[0]
        v2 = coords[2] - coords[0]
        normal = np.cross(v1, v2)
        # Normalize the normal vector
        norm = np.linalg.norm(normal)
        if norm > 1e-6:
            normal = normal / norm
        else:
            # Fallback: use z-axis if ring is degenerate
            normal = np.array([0.0, 0.0, 1.0])
    else:
        normal = np.array([0.0, 0.0, 1.0])

    return centroid, normal


def _is_point_in_cylinder(
    point: np.ndarray,
    ring_centroid: np.ndarray,
    normal: np.ndarray,
    cylinder_radius: float,
    cylinder_height: float,
) -> bool:
    """Check if a point is within a cylinder defined by ring plane and height.

    Args:
        point: 3D coordinates of the point to check.
        ring_centroid: Center of the aromatic ring.
        normal: Normal vector of the ring plane (normalized).
        cylinder_radius: Radius of the cylinder.
        cylinder_height: Height of cylinder extending on either side of plane.

    Returns:
        True if point is within the cylinder, False otherwise.
    """
    # Vector from ring centroid to point
    vec_to_point = point - ring_centroid

    # Project vector onto normal (distance along normal)
    distance_along_normal = np.dot(vec_to_point, normal)

    # Check if point is within cylinder height
    if abs(distance_along_normal) > cylinder_height:
        return False

    # Compute perpendicular distance from axis
    projection = distance_along_normal * normal
    perpendicular_vec = vec_to_point - projection
    perpendicular_distance = np.linalg.norm(perpendicular_vec)

    # Check if within cylinder radius
    return perpendicular_distance <= cylinder_radius


def find_cation_pi_interactions_cylinder(
    points1: List[List[str]],
    points2: List[List[str]],
    cylinder_radius: float,
    cylinder_height: float = 6.0,
) -> List[Tuple[List[str], List[str]]]:
    """Find cation-π interactions using a cylinder-based approach.

    For each aromatic ring (Phe, Tyr, Trp), builds a cylinder extending
    cylinder_height Angstroms on either side of the ring plane. Searches for
    cationic atoms (Lys, Arg, His) within that cylinder.

    Args:
        points1: List of atoms from first chain/structure.
        points2: List of atoms from second chain/structure.
        cylinder_radius: Radius of the cylinder in Angstroms.
        cylinder_height: Height extending on either side of ring plane (default: 6.0).

    Returns:
        List of tuples (aromatic_atom, cation_atom) representing interactions.
    """

    # Group atoms by residue
    def _build_residue_atoms(points: List[List[str]]) -> Dict[Tuple[str, str, str], List[List[str]]]:
        """Group atoms by residue key."""
        by_res: Dict[Tuple[str, str, str], List[List[str]]] = {}
        for atom in points:
            key = _residue_key(atom)
            by_res.setdefault(key, []).append(atom)
        return by_res

    res_atoms1 = _build_residue_atoms(points1)
    res_atoms2 = _build_residue_atoms(points2)

    results = []

    # Process aromatic residues in points1
    for res_key, atoms in tqdm(res_atoms1.items(), leave=False, desc="Processing aromatic rings"):
        res_name = res_key[2].strip().upper()

        # Check if residue is aromatic
        if res_name not in AROMATIC_RING_ATOMS:
            continue

        # Extract ring atoms
        ring_atom_names = AROMATIC_RING_ATOMS[res_name]
        ring_atoms = [atom for atom in atoms if atom[2].strip() in ring_atom_names]

        # Need at least 3 atoms to define a plane
        if len(ring_atoms) < 3:
            continue

        # Compute ring plane (centroid and normal)
        ring_centroid, normal = _compute_ring_plane(ring_atoms)

        # Search for cations in points2
        for res_key2, atoms2 in res_atoms2.items():
            res_name2 = res_key2[2].strip().upper()

            # Check if residue is cationic
            if res_name2 not in CATIONIC_ATOMS:
                continue

            # Extract cationic atoms
            cation_atom_names = CATIONIC_ATOMS[res_name2]
            cation_atoms = [atom for atom in atoms2 if atom[2].strip() in cation_atom_names]

            # Check each cation atom
            for cation_atom in cation_atoms:
                cation_coords = _get_atom_coords(cation_atom)

                # Check if within cylinder
                if _is_point_in_cylinder(
                    cation_coords,
                    ring_centroid,
                    normal,
                    cylinder_radius,
                    cylinder_height,
                ):
                    # Use first ring atom as representative for the interaction
                    # (or could use ring centroid atom representation)
                    results.append((ring_atoms[0], cation_atom))

    # Also process aromatic residues in points2 and cations in points1
    for res_key, atoms in tqdm(res_atoms2.items(), leave=False, desc="Processing aromatic rings (reverse)"):
        res_name = res_key[2].strip().upper()

        if res_name not in AROMATIC_RING_ATOMS:
            continue

        ring_atom_names = AROMATIC_RING_ATOMS[res_name]
        ring_atoms = [atom for atom in atoms if atom[2].strip() in ring_atom_names]

        if len(ring_atoms) < 3:
            continue

        ring_centroid, normal = _compute_ring_plane(ring_atoms)

        # Search for cations in points1
        for res_key2, atoms2 in res_atoms1.items():
            res_name2 = res_key2[2].strip().upper()

            if res_name2 not in CATIONIC_ATOMS:
                continue

            cation_atom_names = CATIONIC_ATOMS[res_name2]
            cation_atoms = [atom for atom in atoms2 if atom[2].strip() in cation_atom_names]

            for cation_atom in cation_atoms:
                cation_coords = _get_atom_coords(cation_atom)

                if _is_point_in_cylinder(
                    cation_coords,
                    ring_centroid,
                    normal,
                    cylinder_radius,
                    cylinder_height,
                ):
                    results.append((ring_atoms[0], cation_atom))

    return results
