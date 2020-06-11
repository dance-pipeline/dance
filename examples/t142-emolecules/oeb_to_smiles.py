"""Converts input OEB file to SMILES file.

Usage:
    python oeb_to_smiles.py INPUT_FILE OUTPUT_FILE
"""

import sys

from openeye import oechem

dataset_stream = oechem.oemolistream(sys.argv[1])
smiles_stream = oechem.oemolostream(sys.argv[2])
for mol in dataset_stream.GetOEMols():
    oechem.OEWriteMolecule(smiles_stream, mol)
smiles_stream.close()
dataset_stream.close()
