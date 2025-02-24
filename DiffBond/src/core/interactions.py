import os
from pathlib import Path
import shutil
import tempfile
import logging

import src.utils.file_handler as file_handler
import src.core.hbondfinder_handler as hbondfinder_handler
import src.utils.distance as distance
from src.utils.graph_handler import convert_to_indexed_edge_list


# Holds code to run the contact mode
def c_mode(PDB_data, dist, verbose=False):
    print("##### Searching contacts within " + str(dist) + "... #####")
    contact_edges = distance.find_electrostatic_interactions(PDB_data[0], PDB_data[1], dist, False)
    atom_index, edge_list = convert_to_indexed_edge_list(contact_edges)
    if verbose:
        print("--------------------------------------------------------------------------------------")
        print("---CONTACT EDGE PREDICTIONS---")
        print("--------------------------------------------------------------------------------------")
        print(edge_list)
    return atom_index, edge_list, contact_edges


# Holds code to run the ionic mode
def i_mode(PDB_data, dist, verbose=False):
    print("##### Searching ionic bonds... #####")
    ionic_edges = distance.find_electrostatic_interactions(PDB_data[0], PDB_data[1], dist, True)
    atom_index, edge_list = convert_to_indexed_edge_list(ionic_edges)
    if verbose:
        print("--------------------------------------------------------------------------------------")
        print("---IONIC BOND PREDICTIONS THAT MEET ELECTROSTATIC CRITERIA---")
        print("--------------------------------------------------------------------------------------")
        print(atom_index)
        print("---Indexed list of edges---")
        print(edge_list)

    return atom_index, edge_list, ionic_edges


def a_mode(PDB_data, dist, bonds, verbose=False):
    print("##### Searching contacts adjacent to bonds #####")
    adj_edges = distance.compareDistAdj(bonds, PDB_data[0], PDB_data[1], dist)
    atom_index, edge_list = convert_to_indexed_edge_list(adj_edges)
    return atom_index, edge_list, adj_edges


def p_mode(PDB_data, dist, ring_radius, verbose=False):
    print("##### Searching cation-pi bonds... #####")
    cation_pi_edges = distance.find_cation_pi_interactions(PDB_data[0], PDB_data[1], dist + 1, ring_radius)
    atom_index, edge_list = convert_to_indexed_edge_list(cation_pi_edges)
    print(cation_pi_edges)


def s_mode(ionic_edges, hbond_edges, verbose=False):
    if ionic_edges is None or hbond_edges is None:
        return

    def match_amino_acid_tuples(list1, list2):
        list1_matches = []
        list2_matches = []
        for pair1 in list1:
            for pair2 in list2:
                paired = False
                # if (pair1[0][5].strip() == pair2[0][5].strip() and pair1[0][4].strip() == pair2[0][4].strip()) and (
                #     pair1[1][5].strip() == pair2[1][5].strip() and pair1[1][4].strip() == pair2[1][4].strip()
                # ):
                #     paired = True
                # if (pair1[0][5].strip() == pair2[1][5].strip() and pair1[0][4].strip() == pair2[1][4].strip()) and (
                #     pair1[1][5].strip() == pair2[0][5].strip() and pair1[1][4].strip() == pair2[0][4].strip()
                # ):
                #     paired = True
                if (pair1[0][5].strip() == pair2[0][5].strip()) and (pair1[1][5].strip() == pair2[1][5].strip()):
                    paired = True
                if (pair1[0][5].strip() == pair2[1][5].strip()) and (pair1[1][5].strip() == pair2[0][5].strip()):
                    paired = True
                # print(pair1, pair2)
                if paired:
                    print(pair1[0][5], pair2[0][5], pair1[1][5], pair2[1][5])
                    if pair1 not in list1_matches:
                        list1_matches.append(pair1)
                    if pair2 not in list2_matches:
                        list2_matches.append(pair2)

        return list1_matches, list2_matches

    # print(ionic_edges, hbond_edges)
    saltbridge_i_edges, saltbridge_h_edges = match_amino_acid_tuples(ionic_edges, hbond_edges)
    atom_index1, edge_list1 = convert_to_indexed_edge_list(saltbridge_i_edges)
    atom_index2, edge_list2 = convert_to_indexed_edge_list(saltbridge_h_edges)

    if verbose:
        print("--------------------------------------------------------------------------------------")
        print("-----------------------------SALT BRIDGE PREDICTIONS-----------------------------------")
        print("--------------------------------------------------------------------------------------")
        # print(saltbridge_i_edges)
        # print(saltbridge_h_edges)
        print(atom_index1)
        print(atom_index2)
        print("---Indexed list of edges---")
        print(edge_list1)
        print(edge_list2)

    return [atom_index1, atom_index2], [edge_list1, edge_list2], [saltbridge_i_edges, saltbridge_h_edges]


# Holds code to run the hbond mode
def h_mode(PDB_data, outputFileName, verbose=False):
    print("##### Searching h-bonds... #####")
    edges_temp = distance.find_electrostatic_interactions(PDB_data[0], PDB_data[1], 3.5, False)

    print("##### Processing h-bonds... #####")

    with tempfile.TemporaryDirectory(delete=False) as temp_dir:
        temp_dir_path = Path(temp_dir)
        HB_file, hb_file = handle_hbond_processing(edges_temp, temp_dir_path, outputFileName)

    if HB_file is None or hb_file is None:
        print("ERROR: No bonds in contact distance to run hbondfinder.")
        return {}, [], []

    # HB_lines = file_handler.parse_file(HB_file, True, 1)
    # hb_lines = file_handler.parse_file(hb_file, True, 1)
    atom_index, edge_list = hbondfinder_handler.parse_hblines_file(hb_file)
    hbond_edges = hbondfinder_handler.edges_with_atom_details(atom_index, edge_list)

    if verbose:
        print("--------------------------------------------------------------------------------------")
        print("---H-BOND PREDICTIONS THAT MEET HBONDFINDER CRITERIA---")
        print("--------------------------------------------------------------------------------------")
        print(atom_index)
        print("---Indexed list of edges---")
        print(edge_list)

    return atom_index, edge_list, hbond_edges


def handle_hbond_processing(edges, temp_dir, output_file_name):
    """
    Process hydrogen bonds and handle related file operations.

    Args:
            edges (list): List of edges found.
            temp_dir (Path): Path to the temporary directory.
            output_file_name (str): The name of the output file.

    Returns:
            tuple: Tuple containing the names of the hydrogen bond files.
    """
    if not edges:
        return None, None

    atoms1 = [edge[0] for edge in edges]
    atoms2 = [edge[1] for edge in edges]

    temp_pdb_path = temp_dir / "temp.pdb"

    # Write atoms to PDB files
    file_handler.write_PDB(str(temp_pdb_path), False, atoms1)
    file_handler.write_PDB(str(temp_pdb_path), True, atoms2)

    output_file_path = Path("results") / output_file_name

    # Run hydrogen bond finder
    try:
        hbondfinder_handler.run_hbondfinder(str(temp_pdb_path))
    except Exception as error:
        logging.error("Failed to run hbondfinder", error)
        shutil.copy(temp_pdb_path, output_file_path / "hbondfinder_error.txt")
        return None, None

    try:
        # Move HBondFinder files
        shutil.move("HBondFinder_temp.txt", output_file_path / "HBondFinder_output.txt")
        shutil.move("hbonds_temp.txt", output_file_path / "hbonds_list.txt")
    except OSError as error:
        logging.error(f"ERROR: {error}")
        return None, None

    return output_file_path / "HbondFinder_output.txt", output_file_path / "hbonds_list.txt"


# TODO: Need to be fixed after refactoring
def handle_interchain_hbond_processing(file):
    """Process interchain hydrogen bonds.

    Args:
            file (str): Path to the PDB file.

    Returns:
            str: Name of the hydrogen bond file.
    """
    try:
        shutil.copyfile(file, "./temp.pdb")
    except OSError:
        logging.error(f"ERROR: Was not able to find {file}")
        return

    if not hbondfinder_handler.run_hbondfinder("temp.pdb"):
        os.remove("temp.pdb")
        return

    hb_file_name = f"{file.split('.')[-2].split('\\')[-1]}.txt"
    try:
        os.rename("HBondFinder_temp.txt", f"HBondFinder{hb_file_name}")
        shutil.move(f"HBondFinder{hb_file_name}", "hbondfinder_data")
        os.rename("hbonds_temp.txt", f"hbonds{hb_file_name}")
        shutil.move(f"hbonds{hb_file_name}", "hbondfinder_data")
    except OSError:
        logging.error("Was not able to move hbond file to hbond_data folder.")
    finally:
        os.remove("temp.pdb")

    return f"HBondFinder{hb_file_name}"
