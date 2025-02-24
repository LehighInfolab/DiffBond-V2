"""
Utility functions for creating and manipulating graphs from PDB data.

Main functions:
- convert_to_indexed_edge_list: Convert PDB tuples to indexed edge list format
- make_graph_hbond: Create graph from hydrogen bond edges
- make_digraph_hbond: Create directed graph from hydrogen bond edges
- visualize_graph: Visualize and save graph plots
"""


import numpy as np
import pandas as pd
import networkx as nx
import pathlib as Path

# from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Any, Optional


def convert_to_indexed_edge_list(pdb_tuples: List[Tuple[List, List]]) -> Tuple[Dict[int, Tuple], List[Tuple[int, int]]]:
    """Convert a list of tuples from a PDB file into an indexed edge list format.

    Args:
        pdb_tuples: Each tuple contains two lists representing bonded atoms.

    Returns:
        Tuple containing:
            - Dictionary mapping indices to atom identifiers
            - List of tuples representing edges between indexed atoms
    """
    atom_index = {}
    edge_list = []
    current_index = 0

    for bond in pdb_tuples:
        for atom in bond:
            atom_id = tuple(atom)  # Use tuple(atom) to create a hashable unique identifier
            if atom_id not in atom_index:
                atom_index[current_index] = atom_id
                current_index += 1

        # Find indices for both atoms in the bond
        atom1_idx = next(idx for idx, atom in atom_index.items() if atom == tuple(bond[0]))
        atom2_idx = next(idx for idx, atom in atom_index.items() if atom == tuple(bond[1]))
        edge_list.append((atom1_idx, atom2_idx))

    return atom_index, edge_list


def write_index_edge_mapping(index: Optional[Dict], edges: Optional[List], name: str, output_dir: Path) -> None:
    """Write node-to-index mapping and edge list to a file.

    Args:
        index: Dictionary mapping node indices to atoms
        edges: List of edges between nodes
        name: Base name for output file
        output_dir: Directory to write output file
    """
    output_file_path = output_dir / (name + "_edges.txt")

    if index is None or edges is None:
        with open(output_file_path, "w") as f:
            f.write("")
        return

    # Convert data to DataFrames for easy CSV writing
    nodes_df = pd.DataFrame(list(index.items()), columns=["Node", "Index"])
    edges_df = pd.DataFrame(edges, columns=["Source", "Target"])

    try:
        with open(output_file_path, "w") as f:
            # f.write("# Node to Index Mapping\n")
            nodes_df.to_csv(f, index=False)
            # f.write("\n# Edges\n")
            edges_df.to_csv(f, index=False)
    except Exception as error:
        raise error


def _process_hbond_node(G: nx.Graph, edge: List, node_idx: int, pos: Dict, color: Dict) -> None:
    """Helper function to process a node in hydrogen bond graphs.

    Args:
        G: NetworkX graph
        edge: Edge data containing node information
        node_idx: Index of the node (0 or 1)
        pos: Dictionary of node positions
        color: Dictionary of node colors
    """
    node_name = f"{edge[node_idx][0]}_{edge[node_idx][1]}_{edge[node_idx][2]}"

    # Add node attributes
    G.nodes[node_name].update({
        "AA": edge[node_idx][2],
        "coord": edge[node_idx][3],
        "chain": edge[node_idx][0]
    })

    # Handle hbond attribute
    hbond_val = G.nodes[node_name].get("hbond")
    if hbond_val is None:
        G.nodes[node_name]["hbond"] = edge[node_idx][4]
    elif hbond_val != edge[node_idx][4]:
        G.nodes[node_name]["hbond"] = ["donor", "acceptor"]

    # Set position and color
    pos[node_name] = np.array([edge[node_idx][3][0], edge[node_idx][3][1]])
    color[node_name] = (0, 1, 1) if node_idx == 0 else (1, 0.3, 0.3)


def make_graph_hbond(edges: List[Any]) -> Tuple[nx.Graph, Dict, Dict]:
    """Create an undirected graph from hydrogen bond edges.

    Args:
        edges: List of hydrogen bond edges

    Returns:
        Tuple containing:
            - NetworkX graph
            - Dictionary of node positions
            - Dictionary of node colors
    """
    G = nx.Graph()
    pos = nx.get_node_attributes(G, "pos")
    color = nx.get_node_attributes(G, "color")

    if not edges:
        return G, pos, color

    for edge in edges:
        node0_name = f"{edge[0][0]}_{edge[0][1]}_{edge[0][2]}"
        node1_name = f"{edge[1][0]}_{edge[1][1]}_{edge[1][2]}"

        G.add_edge(node0_name, node1_name, bond_type=edge[2][0], weight=edge[2][1])

        for i in range(2):
            _process_hbond_node(G, edge, i, pos, color)

    return G, pos, color


def make_digraph_hbond(edges: List[Any]) -> Tuple[nx.DiGraph, Dict, Dict]:
    """Create a directed graph from hydrogen bond edges.

    Args:
        edges: List of hydrogen bond edges

    Returns:
        Tuple containing:
            - NetworkX directed graph
            - Dictionary of node positions
            - Dictionary of node colors
    """
    G = nx.DiGraph()
    pos = nx.get_node_attributes(G, "pos")
    color = nx.get_node_attributes(G, "color")

    if not edges:
        return G, pos, color

    for edge in edges:
        G.add_edge(edge[0][1], edge[1][1], bond_type=edge[2][0], weight=edge[2][1])

        for i in range(2):
            _process_hbond_node(G, edge, i, pos, color)

    return G, pos, color


def make_graph(edges):
    """Util function for making a networkx graph object from a list of edges. Automatically puts into bipartite format.

    Args:
        edges (list): Graph edges formatted into list format

    Returns:
        G: Networkx graph from edges input
        pos: position of nodes
        color: colors for nodes based on which side of the interface node is on
    """
    G = nx.Graph()
    pos = nx.get_node_attributes(G, "pos")
    color = nx.get_node_attributes(G, "color")

    # if edges is empty, return
    if not edges:
        return G, [], []

    # make sure all nodes on left side are from the same chain
    left_chain = edges[0][0][0]
    # left_nodes = []

    for edge in edges:
        same_chain = False
        # if the first node, edge[0], is not on the same chain as left_chain, then switch positions of the nodes.
        # Use flipped boolean to keep track if edges have been switched

        if edge[0][0] != left_chain:
            temp = edge[0]
            edge[0] = edge[1]
            edge[1] = temp

        node0_name = edge[0][0] + "_" + edge[0][1] + "_" + edge[0][2]
        node1_name = edge[1][0] + "_" + edge[1][1] + "_" + edge[1][2]
        # added edge labels here
        G.add_edge(node0_name, node1_name, bond_type=edge[2][0], weight=edge[2][1])

        # adding node labels
        G.nodes[node0_name]["AA"] = edge[0][2]
        G.nodes[node1_name]["AA"] = edge[1][2]

        G.nodes[node0_name]["coord"] = edge[0][3]
        G.nodes[node1_name]["coord"] = edge[1][3]

        G.nodes[node0_name]["chain"] = edge[0][0]
        G.nodes[node1_name]["chain"] = edge[1][0]

        # use atom positions for position in figures
        pos[node0_name] = np.array([edge[0][3][0], edge[0][3][1]])
        pos[node1_name] = np.array([edge[1][3][0], edge[1][3][1]])

        # use edge to indicate color
        if not same_chain:
            color[edge[0][1]] = (0, 1, 1)
            color[edge[1][1]] = (1, 0.3, 0.3)

    # num_nodes = G.number_of_nodes()
    # num_edges = G.number_of_edges()

    # draw position of nodes on one side of interface separate from other side
    # pos = nx.drawing.layout.bipartite_layout(G, left_nodes)

    return G, pos, color


def visualize_graph(G: nx.Graph, pos: Dict, color: Dict, save_dir_name: str) -> None:
    """Visualize and save a graph plot.

    Args:
        G: NetworkX graph to visualize
        pos: Dictionary of node positions
        color: Dictionary of node colors
        save_dir_name: Path to save the plot
    """
    if nx.is_empty(G):
        return

    plt.subplot(121)
    nx.draw(G, with_labels=True, node_color=list(color.values()), pos=pos)

    edge_labels = nx.get_edge_attributes(G, "bond_type")
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)

    plt.savefig(save_dir_name, bbox_inches="tight")
    plt.show()
