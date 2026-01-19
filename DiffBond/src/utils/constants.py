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
    "TRP": {"charge": 0, "type": "aromatic"},
}

# Atoms involved in ionic bonds
CHARGED_ATOMS: Dict[str, List[str]] = {"positive": ["NZ", "NE", "NH1", "NH2"], "negative": ["OD1", "OD2", "OE1", "OE2"]}

# Aromatic ring atoms for each residue type
AROMATIC_RING_ATOMS: Dict[str, List[str]] = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],  # Benzene ring
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],  # Benzene ring (same as PHE)
    "TRP": ["CG", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3"],  # Indole ring (larger ring)
}

# Cationic atoms for each residue type
CATIONIC_ATOMS: Dict[str, List[str]] = {
    "LYS": ["NZ"],
    "ARG": ["NE", "NH1", "NH2"],
    "HIS": ["ND1", "NE2"],  # Both can be protonated
}

# Amino acid single-letter to three-letter code mapping
SINGLE_TO_THREE_LETTER: Dict[str, str] = {
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

PDB_COORDINATE_INDICES = {"x": 6, "y": 7, "z": 8}

# Common molecular measurements and thresholds
MOLECULAR_MEASUREMENTS = {
    "default_max_distance": 5.0,  # Angstroms
    "default_ring_radius": 3.5,  # Angstroms
    # "default_hbond_length": 1.5    # Angstroms
}
