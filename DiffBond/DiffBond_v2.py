"""
DiffBond - A tool for analyzing molecular interactions in protein structures.

Identifies various types of bonds and interactions within PDB structures:
- Ionic bonds
- Hydrogen bonds
- Salt bridges
- Cation-π interactions
- Contact points

Usage:
    diffbond -i [input_files] -m [modes] [-d distance] [-o output]

Arguments:
    -i/--input: One or two PDB files
    -m/--mode: Interaction types to analyze (c,i,h,S,p)
    -d/--distance: Distance threshold (Å)
    -o/--output: Output directory name
"""

from pathlib import Path
import argparse
import shutil
import itertools
from typing import List, Dict, Any

import src.utils.graph_handler as graph_handler
import src.utils.file_handler as file_handler
import src.core.interactions as interactions
from src.utils.constants import MOLECULAR_MEASUREMENTS

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Identify molecular interactions in PDB structures.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-i", "--input",
        nargs="+",
        required=True,
        type=Path,
        help="Input PDB files (1 or 2 files)"
    )

    parser.add_argument(
        "-o", "--output",
        type=Path,
        help="Output directory name"
    )

    parser.add_argument(
        "-m", "--mode",
        nargs="+",
        required=True,
        choices=['c', 'i', 'h', 's', 'p'],
        help="Analysis modes: c=contact, i=ionic, h=hydrogen, S=salt bridge, p=cation-pi"
    )

    parser.add_argument(
        "-d", "--distance",
        type=float,
        default=MOLECULAR_MEASUREMENTS["default_max_distance"],
        help=f"Distance threshold in Å (default: {MOLECULAR_MEASUREMENTS['default_max_distance']})"
    )

    parser.add_argument(
        "-g", "--graph",
        action="store_true",
        help="Output networkx readable files (.adjlist, .gml)"
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print detailed progress information"
    )

    parser.add_argument(
        "-a", "--adjacent",
        choices=['c', 'i', 'h', 's'],
        help="Calculate atoms adjacent to specified bond type"
    )

    return parser.parse_args()

def generate_output_name(input_files: List[Path], distance: float) -> str:
    """Generate default output name from input files."""
    output_name = "Result"

    for file in input_files:
        output_name += f"_{file.stem}"

    if distance != MOLECULAR_MEASUREMENTS["default_max_distance"]:
        output_name += f"_{distance}"

    return output_name

def setup_directories(output_name: str) -> tuple[Path, Path]:
    """Create and return output directories."""
    results_dir = Path("Results") / output_name
    pdb_dir = results_dir / "pdb"

    results_dir.mkdir(parents=True, exist_ok=True)
    pdb_dir.mkdir(parents=True, exist_ok=True)

    return results_dir, pdb_dir

def copy_input_files(input_files: List[Path], pdb_dir: Path) -> None:
    """Copy input PDB files to output directory."""
    for i, file in enumerate(input_files, 1):
        dest = pdb_dir / f"{i}.pdb"
        try:
            shutil.copy(file, dest)
        except Exception as e:
            print(f"Warning: Failed to copy {file} to {dest}: {e}")

def process_pdb_data(input_files: List[Path]) -> List[List[Any]]:
    """Process input PDB files and return data for analysis."""
    if len(input_files) == 1:
        # Split single PDB into chains
        pdb_data = file_handler.parse_PDB_file(input_files[0])
        if not pdb_data:
            raise ValueError(f"No data found in PDB file: {input_files[0]}")

        chains_data = file_handler.split_PDB_chain(pdb_data)
        chains_comb = list(itertools.combinations(chains_data.keys(), 2))

        return [
            [chains_data[c1], chains_data[c2]]
            for c1, c2 in chains_comb
        ]

    else:
        # Process two separate PDB files
        data_list = []
        for file in input_files:
            data = file_handler.parse_PDB_file(file)
            if not data:
                raise ValueError(f"No data found in PDB file: {file}")
            data_list.append(data)
        return [data_list]

def analyze_interactions(
    pdb_data: List[Any],
    modes: List[str],
    distance: float,
    output_name: str,
    verbose: bool
) -> Dict[str, List[Any]]:
    """Analyze molecular interactions based on specified modes."""
    edge_dict = {}

    # Cache for reused calculations
    contact_edges = []
    ionic_edges = []
    hbond_edges = []

    if 'c' in modes:
        atom_index, edge_list, contact_edges = interactions.c_mode(
            pdb_data, distance, verbose
        )
        edge_dict["contact"] = [atom_index, edge_list, contact_edges]

    if 'i' in modes:
        atom_index, edge_list, ionic_edges = interactions.i_mode(
            pdb_data, distance, verbose
        )
        edge_dict["ionic"] = [atom_index, edge_list, ionic_edges]

    if 'p' in modes:
        cation_pi_edges = interactions.p_mode(
            pdb_data,
            distance,
            ring_radius=MOLECULAR_MEASUREMENTS["default_ring_radius"],
            verbose=verbose
        )
        edge_dict["cationpi"] = [atom_index, edge_list, cation_pi_edges]

    if 'h' in modes:
        atom_index, edge_list, hbond_edges = interactions.h_mode(
            pdb_data, output_name, verbose
        )
        edge_dict["hbond"] = [atom_index, edge_list, hbond_edges]

    if 's' in modes:
        # Get ionic and hbond edges if not already calculated
        if not ionic_edges:
            _, _, ionic_edges = interactions.i_mode(pdb_data, distance, verbose)
        if not hbond_edges:
            _, _, hbond_edges = interactions.h_mode(pdb_data, output_name, verbose)

        atom_index, edge_list, saltbridge_edges = interactions.s_mode(
            ionic_edges, hbond_edges, verbose
        )
        edge_dict["saltbridge_ionic"] = [atom_index[0], edge_list[0], saltbridge_edges[0]]
        edge_dict["saltbridge_hbond"] = [atom_index[1], edge_list[1], saltbridge_edges[1]]

    return edge_dict

def main():
    args = parse_arguments()

    # Generate output name if not provided
    output_name = args.output or generate_output_name(args.input, args.distance)

    # Setup directory structure
    results_dir, pdb_dir = setup_directories(output_name)
    print(f"Results will be saved to: {results_dir}")

    # Copy input files
    copy_input_files(args.input, pdb_dir)

    # Process PDB data
    try:
        pdb_data_list = process_pdb_data(args.input)
    except ValueError as e:
        print(f"Error processing PDB data: {e}")
        return

    # Analyze each PDB dataset
    for i, pdb_data in enumerate(pdb_data_list):
        if len(args.input) == 1:
            print(f"\nAnalyzing chain pair {i+1}/{len(pdb_data_list)}")

        # Analyze interactions
        edge_dict = analyze_interactions(
            pdb_data,
            args.mode,
            args.distance,
            output_name,
            args.verbose
        )

        # Process adjacent atoms if requested
        if args.adjacent:
            # Implementation of adjacent atom analysis...
            pass

        # Generate graph files if requested
        if args.graph:
            for edge_type, edge_data in edge_dict.items():
                print(f"Writing graph files for {edge_type}...")
                graph_handler.write_index_edge_mapping(
                    edge_data[0],
                    edge_data[1],
                    edge_type,
                    results_dir
                )
                print(f"Finished writing {edge_type} graph files")


if __name__ == "__main__":
    main()


# python DiffBond_v2.py -i datasets/00001/H1-E/hydrogensAdded.pdb datasets/00001/H2-I/hydrogensAdded.pdb -m h
# python hbondfinder.py -i datasets/00001/H1-E/hydrogensAdded.pdb -j libs/JSON_Files/acceptors_donors_dict.json

# python DiffBond_v2.py -i datasets/1brs_muts/00106/H1-A/final_half1.pdb
# datasets/1brs_muts/00106/H2-D/final_half2.pdb -m i -g
