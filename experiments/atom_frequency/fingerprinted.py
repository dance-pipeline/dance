import matplotlib.pyplot as plt 
import math
import argparse
from collections import defaultdict
from openeye import oechem, oeomega

def makeHistogram():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles_database",
                        type=str,
                        required=True)
    args = parser.parse_args()

    ifs = oechem.oemolistream(args.smiles_database)
    
    histogram = []
    
    for molecule in ifs.GetOEMols():
        mol, status = smiles2oemol(oechem.OEMolToSmiles(molecule))
        length = len([atom for atom in mol.GetAtoms()])
        histogram.append(length)
    plt.hist(histogram, density=True, bins=40)
    plt.ylabel("Frequency")
    plt.xlabel("Atom Count")
    plt.title("Atom Frequency Probability of Fingerprinted Molecules")
    plt.show()

def smiles2oemol(smiles):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)

    omega = oeomega.OEOmega()

    omega.SetMaxConfs(1)
    omega.SetIncludeInput(True)
    omega.SetIncludeInput(True)
    omega.SetSampleHydrogens(True)
    omega.SetStrictStereo(True)
    omega.SetStrictAtomTypes(True)
    omega.SetIncludeInput(False)
    status = omega(mol)
    return (mol, status)

if __name__ == "__main__":
    makeHistogram()