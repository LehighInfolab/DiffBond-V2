"""
Script: DiffBond_v2.py - Main entry point for DiffBond tool for analyzing molecular interactions in protein structures.

Identifies various types of bonds and interactions within PDB structures:
- Ionic bonds
- Hydrogen bonds
- Salt bridges
- Cation-π interactions
- Contact points

Usage:
    Diffbond_v2.py -i [input_files] -m [modes] [-d distance] [-o output]

Arguments:
    -i/--input: One or two PDB files
    -m/--mode: Interaction types to analyze (c,i,h,s,p)
    -d/--distance: Distance threshold (Å)
    -o/--output: Output directory name
"""

from pathlib import Path
import argparse
import itertools
import logging
from typing import List, Dict, Any

import src.utils.logging_handler as logging_handler
import src.utils.graph_handler as graph_handler
import src.utils.file_handler as file_handler
import src.utils.mutation_handler as mutation_handler

import src.core.interactions as interactions
from src.core.interactions import _residue_key
from src.utils.constants import MOLECULAR_MEASUREMENTS
import src.core.hbondfinder_handler as hbondfinder_handler

logging_handler.setup_logging()
logger = logging.getLogger(__name__)


def parse_arguments() -> argparse.Namespace:
    """Parse and validate command line arguments for the tool."""
    parser = argparse.ArgumentParser(
        description="Identify molecular interactions in PDB structures.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Input files: either a single PDB (chains will be split automatically) or two PDBs
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        type=Path,
        help="Input PDB files (1 or 2 files)",
    )

    # Optional custom name for the results directory
    parser.add_argument("-o", "--output", type=Path, help="Output directory name (under ./Results)")

    # Modes: specify which interaction types to analyze
    parser.add_argument(
        "-m",
        "--mode",
        required=True,
        help="Interaction types to analyze, e.g. 'c,i,h,s,p,b' "
        "(c=contact, i=ionic, h=hydrogen, s=salt bridge, p=cation-π [cylinder-based], b=covalent-backbone)",
    )

    # Node granularity: atom-level (default) or residue-level nodes
    parser.add_argument(
        "--node-level",
        choices=["atom", "residue"],
        default="atom",
        help="Choose graph nodes as atoms or amino acid residues (default: atom)",
    )

    # Distance strategy for contact-based residue aggregation
    parser.add_argument(
        "--contact-distance",
        choices=["avg-atom", "centroid"],
        default="avg-atom",
        help=(
            "When using --node-level=residue, compute contact distances as average atom-atom distance "
            "(avg-atom) or distance between residue centroids (centroid)."
        ),
    )

    # Global distance cutoff for geometric interactions
    parser.add_argument(
        "-d",
        "--distance",
        type=float,
        default=MOLECULAR_MEASUREMENTS["default_max_distance"],
        help=f"Distance threshold in Å (default: " f"{MOLECULAR_MEASUREMENTS['default_max_distance']})",
    )

    # Cation-π interaction parameters
    parser.add_argument(
        "--cation-pi-radius",
        type=float,
        default=MOLECULAR_MEASUREMENTS["default_ring_radius"],
        help=(
            f"Cylinder radius for cation-π interactions in Å "
            f"(default: {MOLECULAR_MEASUREMENTS['default_ring_radius']})"
        ),
    )
    parser.add_argument(
        "--cation-pi-height",
        type=float,
        default=6.0,
        help="Cylinder height extending on either side of aromatic ring plane in Å (default: 6.0)",
    )
    parser.add_argument(
        "--no-cylinder",
        action="store_true",
        help="Use distance-based method instead of cylinder-based method for cation-π interactions",
    )

    # Optionally export NetworkX-readable graphs
    parser.add_argument(
        "-g",
        "--graph",
        action="store_true",
        help="Write NetworkX-readable files (.adjlist, .gml)",
    )

    # Verbose logging for debugging and progress tracing
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print detailed progress information",
    )

    # (Adjacent mode removed; reserved for future use.)

    # Integration with automated SKEMPI runner (optional):
    # When provided, we will read the SKEMPI CSV and use --row-index to locate
    # the mutant residue(s) for computing mutant-to-node contact distances.
    parser.add_argument(
        "--skempi-csv",
        type=Path,
        help="Path to SKEMPI CSV (semicolon-delimited) to look up mutation info",
    )
    parser.add_argument(
        "--row-index",
        type=int,
        help="Row index (as used by automated script) into the SKEMPI CSV (header skipped)",
    )

    args = parser.parse_args()

    # Normalize modes (accepts comma-separated string, case-insensitive)
    raw = args.mode.replace(" ", "").lower()
    modes = [m for m in raw.split(",") if m]

    valid = {"c", "i", "h", "s", "p", "b"}
    invalid = set(modes) - valid
    if invalid:
        raise SystemExit(f"Invalid mode(s): {sorted(invalid)}. " f"Valid: {sorted(valid)}")
    args.mode = modes
    return args


def process_pdb_data(input_files: List[Path]) -> List[List[Any]]:
    """Read and structure the PDB data.

    Returns:
        - If one PDB given: a list of chain-pairs (each entry is [chain1, chain2]).
        - If two PDBs given: a single entry [[pdb1, pdb2]].
    """
    if len(input_files) == 1:
        # Single file: split into chains and pairwise combine them
        pdb_data = file_handler.parse_PDB_file(input_files[0])
        if not pdb_data:
            raise ValueError(f"No data found in PDB file: {input_files[0]}")

        chains_data = file_handler.split_PDB_chain(pdb_data)
        chains_comb = list(itertools.combinations(chains_data.keys(), 2))

        return [[chains_data[c1], chains_data[c2]] for c1, c2 in chains_comb]

    else:
        # Two files: load both as separate data
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
    verbose: bool,
    node_level: str,
    contact_distance_strategy: str,
    cation_pi_radius: float = None,
    cation_pi_height: float = None,
    use_cylinder: bool = True,
) -> Dict[str, List[Any]]:
    """Run the requested interaction analyses on the given PDB dataset.

    Args:
        pdb_data: List of two PDB chain/dataset objects.
        modes: List of interaction modes to analyze (e.g., ['c', 'i', 'h', 's', 'p', 'b']).
        distance: Distance threshold in Angstroms.
        output_name: Name for output directory.
        verbose: Enable verbose logging.
        node_level: Granularity level ('atom' or 'residue').
        contact_distance_strategy: Strategy for contact distances ('avg-atom' or 'centroid').
        cation_pi_radius: Cylinder radius for cation-π interactions (default: from constants).
        cation_pi_height: Cylinder height for cation-π interactions (default: 6.0 Å).
        use_cylinder: Whether to use cylinder-based method for cation-π (default: True).

    Returns:
        Dictionary mapping interaction type -> [atom_index, edge_list, raw_edges].
    """
    edge_dict: Dict[str, List[Any]] = {}

    # Cache calculated edges so we don't recompute unnecessarily
    ionic_edges = []
    hbond_edges = []

    if "c" in modes:
        c_idx, c_list, c_edges = interactions.c_mode(
            pdb_data,
            distance,
            verbose,
            node_level=node_level,
            contact_distance_strategy=contact_distance_strategy,
        )
        edge_dict["contact"] = [c_idx, c_list, c_edges]

    # Covalent backbone mode: requires contact nodes to define allowed residues
    if "b" in modes:
        # Ensure contact edges available to derive allowed residue set
        if "contact" not in edge_dict:
            c_idx, _c_list, c_edges = interactions.c_mode(
                pdb_data,
                distance,
                verbose,
                node_level=node_level,
                contact_distance_strategy=contact_distance_strategy,
            )
            edge_dict["contact"] = [c_idx, _c_list, c_edges]

        _, _c_list2, c_raw = edge_dict["contact"]

        # Build allowed residue keys from contact raw edges regardless of node granularity
        allowed_residue_keys = set()

        if c_raw and isinstance(c_raw[0], (list, tuple)):
            first = c_raw[0]
            if len(first) == 2 and isinstance(first[0], (list, tuple)) and len(first[0]) > 6:
                # atom-level pairs
                for a1, a2 in c_raw:
                    allowed_residue_keys.add(_residue_key(a1))
                    allowed_residue_keys.add(_residue_key(a2))
            elif len(first) == 3 and isinstance(first[0], (list, tuple)):
                # residue-level triples (r1, r2, d)
                for r1, r2, _d in c_raw:
                    allowed_residue_keys.add(r1)
                    allowed_residue_keys.add(r2)

        b_idx, b_list, b_edges = interactions.cov_mode(
            pdb_data,
            allowed_residue_keys,
            verbose=verbose,
            node_level="residue",
            # contact_distance_strategy=contact_distance_strategy,
        )
        edge_dict["covalent"] = [b_idx, b_list, b_edges]

    if "i" in modes:
        i_idx, i_list, ionic_edges = interactions.i_mode(
            pdb_data,
            distance,
            verbose,
            node_level=node_level,
        )
        edge_dict["ionic"] = [i_idx, i_list, ionic_edges]

    if "p" in modes:
        # Use provided parameters or defaults
        radius = cation_pi_radius if cation_pi_radius is not None else MOLECULAR_MEASUREMENTS["default_ring_radius"]
        height = cation_pi_height if cation_pi_height is not None else 6.0

        p_idx, p_list, p_edges = interactions.p_mode(
            pdb_data,
            distance,
            ring_radius=radius,
            verbose=verbose,
            use_cylinder=use_cylinder,
            cylinder_height=height,
            node_level=node_level,
        )
        edge_dict["cationpi"] = [p_idx, p_list, p_edges]

    if "h" in modes:
        h_idx, h_list, hbond_edges = interactions.h_mode(
            pdb_data,
            output_name,
            verbose,
            node_level=node_level,
        )
        edge_dict["hbond"] = [h_idx, h_list, hbond_edges]

    if "s" in modes:
        # Salt bridges require both ionic + hbond edges - recompute if missing
        if not ionic_edges:
            _, _, ionic_edges = interactions.i_mode(pdb_data, distance, verbose, node_level=node_level)
        if not hbond_edges:
            _, _, hbond_edges = interactions.h_mode(pdb_data, output_name, verbose, node_level=node_level)

        s_idx_pair, s_list_pair, s_edges_pair = interactions.s_mode(ionic_edges, hbond_edges, verbose)
        edge_dict["saltbridge_ionic"] = [s_idx_pair[0], s_list_pair[0], s_edges_pair[0]]
        edge_dict["saltbridge_hbond"] = [s_idx_pair[1], s_list_pair[1], s_edges_pair[1]]

    # Adjacent mode removed

    return edge_dict


def main():
    args = parse_arguments()

    # Determine output folder name (user-specified or generated)
    output_name = args.output or file_handler.generate_output_name(
        args.input, args.distance, MOLECULAR_MEASUREMENTS["default_max_distance"]
    )

    # Best-effort cleanup of any stray HBondFinder temp outputs
    try:
        hbondfinder_handler.cleanup_temp_outputs()
    except Exception as error:
        logger.warning(error)

    results_dir, pdb_dir = file_handler.setup_results_directories(output_name)
    logger.info(f"Results will be saved to: {results_dir}")

    # Save copies of input files to results directory
    file_handler.copy_input_files(args.input, pdb_dir)

    # Parse PDB files into in-memory structures
    try:
        pdb_data_list = process_pdb_data(args.input)
    except ValueError as e:
        logger.error("Error processing PDB data: %s", e)
        return

    # Analyze each dataset (each chain-pair or PDB-pair)
    for i, pdb_data in enumerate(pdb_data_list):
        if len(args.input) == 1:
            logger.info(f"Analyzing chain pair {i + 1} {len(pdb_data_list)}")
        else:
            logger.info(f"Analyzing {args.input}")

        edge_dict = analyze_interactions(
            pdb_data,
            args.mode,
            args.distance,
            output_name,
            args.verbose,
            node_level=args.node_level,
            contact_distance_strategy=args.contact_distance,
            cation_pi_radius=args.cation_pi_radius,
            cation_pi_height=args.cation_pi_height,
            use_cylinder=not args.no_cylinder,
        )

        # TODO: Implement adjacency analysis
        # if args.adjacent: ...

        # Optionally export graph representations
        if args.graph:
            logger.info("Writing graph files...")
            graph_handler.write_graphs_from_edge_dict(edge_dict, results_dir)
            logger.info("Finished writing graph files")

        # Write intramolecular contacts constrained to nodes present in contact edges
        if "contact" in edge_dict:
            try:
                c_idx, _c_edges_with_d, c_raw = edge_dict["contact"]
                out_dir = Path("Results") / output_name
                intra_edges = interactions.compute_intra_contact_edges(
                    pdb_data=pdb_data,
                    contact_index=c_idx,
                    contact_raw=c_raw,
                    max_distance=args.distance,
                    node_level=args.node_level,
                    contact_distance_strategy=args.contact_distance,
                )
                out_file = graph_handler.write_intra_contact_edges(
                    intra_edges, c_idx, out_dir, "intracontact_edges.txt"
                )
                if out_file is not None:
                    logger.info("Wrote intramolecular contact edges: %s", out_file)
            except Exception as _e:
                logger.exception(
                    "Failed computing intramolecular contacts within contact nodes: %s",
                    _e,
                )

        # If SKEMPI CSV info is provided, compute mutant->ALL residue nodes distances (residue level)
        if args.skempi_csv and args.row_index is not None:
            try:
                # Recompute contact edges at residue level to get residue keys
                _, _, res_residue_edges = interactions.c_mode(
                    pdb_data,
                    args.distance,
                    verbose=False,
                    node_level="residue",
                    contact_distance_strategy=args.contact_distance,
                )

                # Read mutation specifications from SKEMPI CSV
                mut_specs = mutation_handler.read_skempi_mutations(args.skempi_csv, args.row_index)

                if not mut_specs:
                    logger.warning(
                        f"No mutation found in SKEMPI CSV at row_index={args.row_index}. "
                        f"CSV file: {args.skempi_csv}. "
                        f"Check if column 4 (index 3) contains a valid mutation token at that row."
                    )
                else:
                    # Collect all residue keys from interaction outputs
                    union_residue_keys = mutation_handler.collect_residue_keys_from_edge_dict(
                        edge_dict, res_residue_edges
                    )

                    if not union_residue_keys:
                        logger.warning(
                            f"No residue keys found in interaction outputs. "
                            f"Cannot compute mutant distances. "
                            f"Interactions computed: {list(edge_dict.keys())}"
                        )
                    else:
                        # Compute distances from mutants to all residue nodes
                        mut_index_map, mut_edges = mutation_handler.compute_mutant_distances(
                            pdb_data,
                            mut_specs,
                            union_residue_keys,
                            args.contact_distance,
                        )

                        # Write indexed edge graph in same format as other interactions
                        if mut_edges:
                            out_dir = Path("Results") / output_name
                            out_dir.mkdir(parents=True, exist_ok=True)
                            try:
                                graph_handler.write_index_edge_mapping(
                                    index=mut_index_map,
                                    edges=mut_edges,
                                    name="mutant",
                                    output_dir=out_dir,
                                )
                                # Count unique mutation nodes (they have indices 0, 1, 2, ...)
                                num_mutations = len(set(src for src, _, _ in mut_edges))
                                logger.info(
                                    f"Wrote mutant edge graph: {len(mut_edges)} edges "
                                    f"from {num_mutations} mutant node(s) to {len(union_residue_keys)} residue nodes"
                                )
                            except Exception as e:
                                logger.exception("Failed writing mutant edge graph: %s", e)
                        else:
                            logger.warning(
                                f"No mutant edges computed. Mutation token(s): {mut_specs}. "
                                f"Possible reasons: "
                                f"(1) Mutation token format invalid (expected format like 'LI45G'), "
                                f"(2) Mutant residue not found in PDB structure, "
                                f"(3) No atoms found for mutant residue, or "
                                f"(4) All target residues were filtered out."
                            )
            except Exception as err:
                logger.exception("Mutant distance computation failed: %s", err)


if __name__ == "__main__":
    main()

# Example commands:
# python DiffBond_v2.py -i datasets/1brs_muts/00106/H1-A/final_half1.pdb \
#   datasets/1brs_muts/00106/H2-D/final_half2.pdb -m h
# python src/scripts/hbondfinder.py -i datasets/1brs_muts/00106/H1-A/final_half1.pdb \
#   -j src/JSON_Files/acceptors_donors_dict.json
# python DiffBond_v2.py -i datasets/1brs_muts/00106/H1-A/final_half1.pdb \
#   datasets/1brs_muts/00106/H2-D/final_half2.pdb -m i -g
