"""Analyze the eMolecules t142 dataset.

Saves CDF plots of each fingerprint component to `fingerprint_cdf.png`

Usage:
    python dataset_analysis.py --dataset-oeb [FILE]
"""
import argparse
from pprint import pprint

import matplotlib.pyplot as plt
from openeye import oechem

from dance import DancePipeline


def main():
    """Parses arguments and prints out the analysis."""
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--dataset-oeb", default="dataset.oeb")
    args = parser.parse_args()

    # Retrieve the fingerprints from the molecules.
    dataset_stream = oechem.oemolistream(args.dataset_file)
    fingerprints = []
    for mol in dataset_stream.GetOEMols():
        fingerprints.append(DancePipeline.get_fingerprint_from_mol(mol))
    dataset_stream.close()

    # Plot CDFs of the fingerprints, display them, and save them to a file.
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    num_atoms = sorted(i[0] for i in fingerprints if i[1] > -1)
    pprint(num_atoms)
    ax[0].set_title("CDF of Molecule Size (Number of atoms)")
    ax[0].hist(num_atoms, bins=num_atoms, cumulative=True)
    ax[0].set_xlabel("Number of atoms")
    ax[0].set_ylabel("Molecules")
    if len(fingerprints) == 20: ax[0].set_yticks(range(0, 21, 2))
    ax[0].grid(True)

    central_wbo = sorted(i[1] for i in fingerprints if i[1] > -1.0)
    pprint(central_wbo)
    ax[1].set_title("CDF of Central WBO")
    ax[1].hist(central_wbo, bins=central_wbo, cumulative=True)
    ax[1].set_xlabel("Central WBO")
    ax[1].set_ylabel("Molecules")
    if len(fingerprints) == 20: ax[0].set_yticks(range(0, 21, 2))
    ax[1].grid(True)

    fig.savefig("results/fingerprint_cdf.png")
    plt.show()


if __name__ == "__main__":
    main()
