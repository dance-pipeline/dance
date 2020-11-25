import matplotlib.pyplot as plt 
import math
import argparse
from collections import defaultdict
from openeye import oechem

def makeHistogram():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles_database",
                        type=str,
                        required=True)
    args = parser.parse_args()

    ifs = oechem.oemolistream(args.smiles_database)
    
    histogram = []
    
    for molecule in ifs.GetOEMols():
        smiles_mol = oechem.OEMolToSmiles(molecule)
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smiles_mol)
        length = len([atom for atom in mol.GetAtoms()])
        histogram.append(length)
    plt.hist(histogram, color="red", density=True, bins=40)
    plt.ylabel("Frequency")
    plt.xlabel("Atom Count")
    plt.title("Atom Frequency Probability of Randomly Selected Molecules")
    plt.show()

if __name__ == "__main__":
    makeHistogram()
