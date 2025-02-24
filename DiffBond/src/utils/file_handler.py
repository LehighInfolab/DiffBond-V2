"""
Utility functions for parsing and writing PDB files and handling file operations.

Main functions:
- parse_PDB_file: Parse a PDB file into a structured format
- write_PDB: Write atom data to PDB format
- make_output_dir: Create output directory for results
"""

import os
from pathlib import Path
from typing import List, Dict, Union


def parse_PDB_file(file: Union[str, Path]) -> List[List[str]]:
    """Parse a PDB file and extract atom information.

    Args:
        file: Path to the PDB file

    Returns:
        List of parsed atom data where each atom is represented as a list of properties
    """
    try:
        with open(file, "r") as f:
            # Only process ATOM lines
            atom_lines = [line for line in f if line.startswith("ATOM")]

        all_lines = []
        for line in atom_lines:
            parsed_line = [
            line[0:6].strip(),  # Record name
            line[6:11].strip(),  # Atom serial number
            line[12:16].strip(),  # Atom name
            # line[16].strip(),  # Alternate location indicator
            line[17:20].strip(),  # Residue name
            line[21].strip(),  # Chain identifier
            line[22:26].strip(),  # Residue sequence number
            # line[26].strip(),  # Code for insertion of residues
            line[30:38].strip(),  # x coordinate
            line[38:46].strip(),  # y coordinate
            line[46:54].strip(),  # z coordinate
            line[54:60].strip(),  # Occupancy
            line[60:66].strip(),  # Temperature factor
            line[76:78].strip(),  # Element symbol
            # line[78:80].strip(),  # Charge
        ]

            # Standardize 4-character residue names
            if len(parsed_line[3]) == 4:
                parsed_line[3] = parsed_line[3][1:]

            all_lines.append(parsed_line)

        return all_lines

    except Exception as e:
        raise RuntimeError(f"Error parsing PDB file {file}: {str(e)}")




def write_PDB(output_name: str, append: bool, atoms_list: List[List[str]]) -> None:
    """Write atom data to a PDB file.

    Args:
        output_name: Path to output file
        append: If True, append to existing file; if False, create new file
        atoms_list: List of atom data to write
    """
    mode = "a" if append else "w"
    with open(output_name, mode) as f:
        for atom in atoms_list:
            # Format each field according to PDB specifications
            record = atom[0].ljust(4)
            atom_num = atom[1].rjust(5)
            atom_name = atom[2].rjust(4)
            res_name = atom[3][1:].ljust(3) if len(atom[3]) == 4 else atom[3].ljust(3)
            chain = atom[4].rjust(1)
            res_num = atom[5].rjust(4)
            x = f"{float(atom[6]):8.3f}".rjust(8)
            y = f"{float(atom[7]):8.3f}".rjust(8)
            z = f"{float(atom[8]):8.3f}".rjust(8)

            line = f"{record}  {atom_num} {atom_name} {res_name} {chain}{res_num}    {x}{y}{z}\n"
            f.write(line)


def split_PDB_chain(PDB_data: List[List[str]]) -> Dict[str, List[List[str]]]:
    """Split PDB data by chain identifier.

    Args:
        PDB_data: List of parsed atom data

    Returns:
        Dictionary mapping chain IDs to their atom data
    """
    chains: Dict[str, List[List[str]]] = {}
    for line in PDB_data:
        chain = str(line[4])
        if chain not in chains:
            chains[chain] = []
        chains[chain].append(line)
    print("All chains in PDB:", chains.keys())
    return chains


def make_output_dir(output_path: Path) -> Path:
    """Create output directory if it doesn't exist.

    Args:
        output_path: Path to create

    Returns:
        Created Path object
    """
    output_path.mkdir(parents=True, exist_ok=True)
    return output_path


def make_chain_comb_dir(output_file_name: str, chains: List[str]) -> str:
    """Create directory for chain combination results.

    Args:
        output_file_name: Base name for output directory
        chains: List of chain identifiers

    Returns:
        Path to created directory
    """
    chain_dir = Path("Results") / output_file_name / f"{chains[0]}_{chains[1]}"
    chain_dir.mkdir(parents=True, exist_ok=True)
    return str(chain_dir)
