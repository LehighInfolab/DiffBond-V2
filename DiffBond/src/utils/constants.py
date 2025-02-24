"""
Constants used across the DiffBond package.

Contains definitions for:
- Amino acid properties and classifications
- Atom types and properties
- Common molecular measurements
"""

from typing import Dict, List, Any

# Amino acid classifications and properties
AMINO_ACID_PROPERTIES: Dict[str, Dict[str, Any]] = {
    # Charged residues with their charge values
    "ARG": {"charge": 1, "type": "cationic"},
    "LYS": {"charge": 1, "type": "cationic"},
    "HIS": {"charge": 1, "type": None},
    "ASP": {"charge": -1, "type": None},
    "GLU": {"charge": -1, "type": None},

    # Aromatic residues
    "PHE": {"charge": 0, "type": "aromatic"},
    "TYR": {"charge": 0, "type": "aromatic"},
    "TRP": {"charge": 0, "type": "aromatic"}
}

# Atoms involved in ionic bonds
CHARGED_ATOMS: Dict[str, List[str]] = {
    "positive": ["NZ", "NE", "NH1", "NH2"],
    "negative": ["OD1", "OD2", "OE1", "OE2"]
}

PDB_COORDINATE_INDICES = {
    "x": 6,
    "y": 7,
    "z": 8
}

# Common molecular measurements and thresholds
MOLECULAR_MEASUREMENTS = {
    "default_max_distance": 5.0,  # Angstroms
    # "default_ring_radius": 3.5,   # Angstroms
    # "default_hbond_length": 1.5    # Angstroms
}