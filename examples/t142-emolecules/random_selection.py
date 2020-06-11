"""Randomly selects N molecules from an input file.

Usage:
    python random_selection.py -n [N] --input [FILE] --output [FILE]
"""

import argparse
import random

from openeye import oechem


def main():
    """Reads in the file and does the selection."""
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-n",
        type=int,
        default=20,
        help="Number of molecules to select.",
    )
    parser.add_argument(
        "--input",
        type=str,
        help="Input molecule file.",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output molecule file.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed.",
    )
    args = parser.parse_args()

    random.seed(args.seed)

    # Input the molecules.
    molecules = []
    input_stream = oechem.oemolistream(args.input)
    for mol in input_stream.GetOEMols():
        molecules.append(oechem.OEMol(mol))

    # Choose the molecules.
    selected_molecules = random.sample(molecules, args.n)

    # Output the molecules.
    output_stream = oechem.oemolostream(args.output)
    for mol in selected_molecules:
        oechem.OEWriteMolecule(output_stream, mol)


if __name__ == "__main__":
    main()
