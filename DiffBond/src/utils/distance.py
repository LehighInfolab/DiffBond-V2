"""
Utility functions for calculating distances and finding molecular interactions.

Main functions:
- find_electrostatic_interactions: Find ionic bonds between charged residues
- find_cation_pi_interactions: Find interactions between aromatic rings and cations
- find_nearby_contacts: Find atoms within specified distance of contact points
"""

from typing import List, Tuple
import math
import numpy as np
from tqdm import tqdm

from src.utils.constants import (
    AMINO_ACID_PROPERTIES,
    CHARGED_ATOMS,
    PDB_COORDINATE_INDICES as PDB_IDX
)

def euclidean_distance(point1: List[str], point2: List[str]) -> float:
    """Calculate 3D Euclidean distance between two points.

    Expects points in PDB format where coordinates are at indices defined in PDB_COORDINATE_INDICES.
    """
    return math.sqrt(sum(
        (float(point1[PDB_IDX[coord]]) - float(point2[PDB_IDX[coord]])) ** 2
        for coord in ["x", "y", "z"]
    ))

def find_nearby_contacts(
    contact_points: List[Tuple[List[str], List[str]]],
    points1: List[List[str]],
    points2: List[List[str]],
    max_distance: float
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
    check_charge: bool = False
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
    ring_radius: float
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