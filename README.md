# DiffBond

## Purpose:

DiffBond identifies and classifies three intermolecular bonds in protein complexes: ionic bonds, salt bridges, and hydrogen bonds.

Input: A protein complex separated into two PDB files for each half of the complex.
Output: A set of files representing the bond network at the interface.
    - A folder copying the pdb files for downstream analysis
    - A .gml file: JSON file of network with 1. Amino acid as nodes, 2. Two amino acids that form a bond as edges, and related pdb features.
    - A .adjlist file: Adjlist format of graph network

    - Files can be read back into a networkx graph object for ease of use

## Command-Line Options:

DiffBond
	- This file contains code for most of the main search functions for finding atoms that meet criteria to form ionic bonds, hydrogen bonds, and salt bridges within a PDB structure.
	
	Usage: -i [ Input files ] -m [ Search Mode ] -d [ Distance threshold ] -o [ Output file name]
	
	-i          Can take 1 or 2 input files. If 1 input provided file, DiffBond will split PDB into chains and find inter-chain bonds. If 2 input provided, DiffBond will find intermolecular bonds between the two input files.

	-m        Search mode. Can be multiple combinations of the following options: Contact = c, Ionic bond = i, Hydrogen bond = h, Salt bridge = S, Cation pi = p. Must include at least 1 option.

	-d         Distance threshold for search distances between atoms, in angstrom units. 

	-o         Output file name.

## Summary of executables and lib files

- DiffBond_v2.py : Executable to generate graph network files for one protein complex

- DiffBond.py : Old version of DiffBond

- hbondfinder.py : Executable for calculating hydrogen bonds - required in the same folder of DiffBond_v2.py in order to calculate hydrogen bonds

- reformatPDBtoSKA.py : Short script to reformat PDB files to be compatible with ska (ska requires all amino acids to be 'A' chain)

- run_diffbond.py : Script to execute multiple instances of DiffBond_v2.py on the SKEMPI dataset or a dataset folder.

- TargetSearchV2.py : DEPRECATED - Helper file for DiffBond_v2.py to search PDB files for bonds.

- lib/
    - graph_utils.py : Utility functions to parse and work with networkx graphs.

    - hbondfinder_utils.py : Utility functions to parse hbondfinder outputs into a list containing hbond edges

    - PDB_HB_parser.py : Miscellaneous utility functions for PDB parsing

    - PDBGreedySearch.py : DEPRECATED

- JSON_Files/
    - Two JSON files mapping amino acid 3-letter codes to its corresponding atom identifier.

- skempi_toy : Toy dataset with a few skempi PDB complexes for testing.