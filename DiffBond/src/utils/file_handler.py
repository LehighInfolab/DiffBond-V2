"""
Utility functions for parsing and writing PDB files and handling file operations.

Main functions:
- parse_PDB_file: Parse a PDB file into a structured format
- write_PDB: Write atom data to PDB format
- make_output_dir: Create output directory for results
"""

from pathlib import Path
from typing import List, Dict, Union
import shutil


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
            # Extract residue sequence number (columns 22-25) and insertion code (column 26)
            resseq_num = line[22:26].strip()  # Residue sequence number (right-justified, 4 chars)
            insertion_code = line[26:27].strip() if len(line) > 26 else ""  # Insertion code (column 26)

            # Combine residue number with insertion code if present (e.g., "60" + "B" -> "60B")
            resseq = resseq_num + insertion_code.lower() if insertion_code else resseq_num

            parsed_line = [
                line[0:6].strip(),  # Record name
                line[6:11].strip(),  # Atom serial number
                line[12:16].strip(),  # Atom name
                # line[16].strip(),  # Alternate location indicator
                line[17:20].strip(),  # Residue name
                line[21].strip(),  # Chain identifier
                resseq,  # Residue sequence number with insertion code (e.g., "60b")
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
            # PDB ATOM line format (columns):
            # 0-5: Record name, 6-10: Atom serial, 12-15: Atom name, 17-19: Residue name
            # 21: Chain, 22-25: Residue number, 26: Insertion code, 30-37: x, 38-45: y, 46-53: z

            record = atom[0].ljust(6)  # "ATOM  " (6 chars)
            atom_num = atom[1].rjust(5)  # Right-justified in columns 6-10
            # Left-justified in columns 17-19
            res_name = atom[3][1:].ljust(3) if len(atom[3]) == 4 else atom[3].ljust(3)
            chain = atom[4] if atom[4] else " "  # Column 21 (1 char, space if empty)

            # Handle residue number and insertion code
            # atom[5] may contain combined residue number + insertion code (e.g., "60B", "100a")
            resseq_str = str(atom[5]).strip()

            # Check if last character is a letter (insertion code)
            if resseq_str and resseq_str[-1].isalpha():
                # Split residue number and insertion code
                # Residue number (columns 22-25, right-justified)
                res_num = resseq_str[:-1].rjust(4)
                insertion_code = resseq_str[-1].upper()  # Insertion code (column 26, 1 char)
            else:
                # No insertion code - use space in column 26
                # Residue number (columns 22-25, right-justified)
                res_num = resseq_str.rjust(4)
                insertion_code = " "  # Space in column 26 when no insertion code

            # Format coordinates with 8.3f format (8 chars total, 3 decimal places)
            x = f"{float(atom[6]):8.3f}"  # Columns 30-37
            y = f"{float(atom[7]):8.3f}"  # Columns 38-45
            z = f"{float(atom[8]):8.3f}"  # Columns 46-53

            # Build PDB line with exact column positions
            # PDB format: ATOM  serial name alt resname chain resseq icode    x      y      z
            # Columns:    0-5   6-10  12-15 16  17-19   21   22-25  26   27-29 30-37 38-45 46-53
            # Note: atom_name should be left-justified in columns 12-15, not right-justified
            atom_name_formatted = atom[2].ljust(4)  # Left-justified in columns 12-15
            line = (
                f"{record}"  # 0-5: "ATOM  "
                f"{atom_num}"  # 6-10: Atom serial (5 chars)
                f" "  # 11: Space
                f"{atom_name_formatted}"  # 12-15: Atom name (4 chars, left-justified)
                f" "  # 16: Alternate location indicator (space)
                f"{res_name}"  # 17-19: Residue name (3 chars)
                f" "  # 20: Space
                f"{chain}"  # 21: Chain identifier (1 char)
                f"{res_num}"  # 22-25: Residue number (4 chars, right-justified)
                f"{insertion_code}"  # 26: Insertion code (1 char)
                f"   "  # 27-29: Three spaces
                f"{x}"  # 30-37: x coordinate (8 chars)
                f"{y}"  # 38-45: y coordinate (8 chars)
                f"{z}"  # 46-53: z coordinate (8 chars)
                f"\n"
            )
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


def pretty_distance(d: float) -> str:
    """Convert a float distance to a compact string form for filenames.

    Examples:
        2.5 -> "2p5"
        3.0 -> "3"

    Args:
        d: Distance value

    Returns:
        Compact string representation
    """
    s = f"{d:.2f}".rstrip("0").rstrip(".")
    return s.replace(".", "p")


def generate_output_name(input_files: List[Path], distance: float, default_distance: float) -> str:
    """Generate a default results folder name from input files and distance.

    Args:
        input_files: List of input PDB file paths
        distance: Distance threshold used
        default_distance: Default distance threshold to check against

    Returns:
        Output directory name
    """
    parts = ["Result"]
    parts.extend(f"{f.stem}" for f in input_files)

    # Append non-default distance to differentiate runs
    if distance != default_distance:
        parts.append(pretty_distance(distance) + "A")
    return "_".join(parts)


def setup_results_directories(output_name: str) -> tuple[Path, Path]:
    """Create Results/<output_name>/ and subdirectories for PDBs.

    Args:
        output_name: Name for the output directory

    Returns:
        Tuple of (results_dir, pdb_dir)
    """
    results_dir = Path("Results") / output_name
    pdb_dir = results_dir / "pdb"

    results_dir.mkdir(parents=True, exist_ok=True)
    pdb_dir.mkdir(parents=True, exist_ok=True)

    return results_dir, pdb_dir


def copy_input_files(input_files: List[Path], pdb_dir: Path) -> None:
    """Copy the original input PDBs into the results folder for recordkeeping.

    Args:
        input_files: List of input PDB file paths
        pdb_dir: Destination directory for copied files
    """
    for i, file in enumerate(input_files, 1):
        dest = pdb_dir / f"{i}.pdb"
        try:
            shutil.copy(file, dest)
        except Exception as e:
            print(f"Warning: Failed to copy {file} to {dest}: {e}")
