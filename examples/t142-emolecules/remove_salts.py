"""Removes salts from the molecules in a given input SMILES file.

The salts are indicated by sections of the SMILES string after a "."

Usage:
    python remove_salts.py [INPUT_SMILES] [OUTPUT_SMILES]
"""
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python remove_salts.py [INPUT_SMILES] [OUTPUT_SMILES]")
        sys.exit(1)

    with open(sys.argv[1], "r") as original_molecules:
        with open(sys.argv[2], "w") as molecules_no_salt:
            for line in original_molecules:
                line = line.strip()

                # Only take the first molecule, not the salt
                line = line.split(".")[0]

                molecules_no_salt.write(line + "\n")
