# DiffBond

## Purpose:

DiffBond identifies and classifies three intermolecular bonds in protein complexes: ionic bonds, salt bridges, and hydrogen bonds.

## Command-Line Options:

-i \
 Input PDB files to be compared if multiple. This option will default to finding intermolecular bonds at the interface between input PDBs.

-d \
Resolution for distance checking. Increasing distance will lessen strictness in search, and vice versa. Default is 5 angstroms.

-m \
Search mode option. Allows you to perform different searches instead of the default (all three of ionic bond, hydrogen bond, salt bridge search). Contact search = [], ionic bond search = i, hydrogen bond search =h, salt bridge search = s.

-o \
Name of output file. Default is Contact*[PDB1]*[PDB2].pdb.
