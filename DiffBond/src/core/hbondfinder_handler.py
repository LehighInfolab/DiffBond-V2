"""
Utility functions for working with HBondFinder program and its output.

Main functions:
- run_hbondfinder: Execute HBondFinder on a PDB file
- parse_hbond_lines: Parse HBondFinder output into edge format
- parse_hblines_file: Parse HBondFinder file into atom and edge data
"""

import os
import logging
import subprocess
from pathlib import Path
from typing import List, Dict, Tuple, Optional

logger = logging.getLogger(__name__)

"""
This file contains some utility functions for working with hbondfinder format
    - run hbondfinder in one function
    - make a folder for holding hbondfinder data
    - batch run hbondfinder (not yet implemented)
    - provide an input of HBondFinder file lines split using spaces,
    and will parse it to find intramolecular or intermolecular edges between atoms
"""


def get_project_root() -> Path:
    """Get the project root directory.
    Assumes this file is in DiffBond/src/core/
    """
    return Path(__file__).parent.parent.parent


def run_hbondfinder(file: Path | str) -> bool:
    """Run HBondFinder on an input PDB file.

    Uses acceptors_donors_dict.json as mapping for amino acids to acceptor/donor atoms.

    Args:
        file: Path to PDB file to analyze

    Returns:
        True if HBondFinder ran successfully, False if input file was empty

    Raises:
        subprocess.CalledProcessError: If HBondFinder execution fails
        FileNotFoundError: If required files are not found
    """
    file = Path(file) if isinstance(file, str) else file

    if not file.exists():
        raise FileNotFoundError(f"Input PDB file not found: {file}")
    if file.stat().st_size == 0:
        return False

    root_dir = get_project_root()

    hbondfinder_script = root_dir / "src" / "scripts" / "hbondfinder.py"
    json_file = root_dir / "src" / "JSON_Files" / "acceptors_donors_dict.json"

    if not hbondfinder_script.exists():
        raise FileNotFoundError(f"HBondFinder script not found at {hbondfinder_script}")
    if not json_file.exists():
        raise FileNotFoundError(f"Acceptors/donors JSON file not found at {json_file}")

    command = ["python", str(hbondfinder_script), "-i", str(file), "-j", str(json_file)]
    logger.info(f"Running hbondfinder.py on {file.name} using command ' {command} ' ...")
    subprocess.run(command, check=True)
    return True


def parse_hbond_lines(lines: List[List[str]], intermolecular: bool = False) -> List[List[str]]:
    """Parse HBondFinder output lines into edge format.

    Args:
        lines: Lines from HBondFinder output file split by spaces
        intermolecular: If True, only include bonds between different chains

    Returns:
        List of edges where each edge is [donor, acceptor, distance]
    """
    edges = []
    for line in lines:
        same_chain = line[0] == line[4]

        if not intermolecular or not same_chain:
            donor = line[0] + line[1]
            acceptor = line[4] + line[5]
            dist = line[8]
            edges.append([donor, acceptor, dist])

    return edges


def cleanup_temp_outputs() -> None:
    """HBondFinder temp output files in CWD."""
    try:
        for fname in ("HBondFinder_temp.txt", "hbonds_temp.txt"):
            if os.path.isfile(fname):
                Path(fname).unlink()
    except Exception:
        # Silently ignore cleanup errors
        pass


def parse_hblines_file(
    file_path: Path,
) -> Tuple[Dict[int, List[str]], List[Tuple[int, int]]]:
    """Parse a HBondFinder output file into atom and edge data.

    Args:
        file_path: Path to HBondFinder output file

    Returns:
        Tuple containing:
            - Dictionary mapping atom indices to atom information
            - List of edges between atoms as (source, target) index pairs

    Raises:
        RuntimeError: If there's an error parsing the file
        FileNotFoundError: If the file doesn't exist
    """
    if not file_path.exists():
        raise FileNotFoundError(f"HBondFinder output file not found: {file_path}")

    try:
        with open(file_path, "r") as file:
            lines = file.readlines()

            # Parse number of atoms
            num_atoms = int(lines[0].split()[1])

            # Parse atom information
            atom_index = {}
            for i in range(1, num_atoms + 1):
                fields = lines[i].split()
                atom_info = [
                    "ATOM",  # Record type
                    fields[4],  # Atom serial number
                    fields[3],  # Atom name
                    fields[1],  # Residue name
                    fields[0],  # Chain ID
                    fields[2],  # Residue number
                    fields[5],  # x coordinate
                    fields[6],  # y coordinate
                    fields[7],  # z coordinate
                    "1.00",  # Occupancy
                    "0.00",  # Temperature factor
                ]
                atom_index[int(fields[4])] = atom_info

            # Parse hydrogen bonds
            num_hbonds = int(lines[num_atoms + 1].split()[1])
            edge_list = []
            for i in range(num_atoms + 2, num_atoms + 2 + num_hbonds):
                source, target = map(int, lines[i].split()[:2])
                edge_list.append((source, target))

            return atom_index, edge_list

    except Exception as e:
        raise RuntimeError(f"Error parsing HBondFinder file {file_path}: {str(e)}")


def edges_with_atom_details(
    atom_index: Dict[int, List[str]], edge_list: List[Tuple[int, int]]
) -> List[Tuple[List[str], List[str]]]:
    """Convert edge list to detailed format including atom information.

    Args:
        atom_index: Dictionary mapping indices to atom information
        edge_list: List of edges as (source, target) index pairs

    Returns:
        List of edges with full atom details for both atoms
    """
    return [(atom_index[edge[0]], atom_index[edge[1]]) for edge in edge_list]
