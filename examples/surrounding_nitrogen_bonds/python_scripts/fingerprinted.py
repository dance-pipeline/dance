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
    ofs = oechem.oemolostream("fingerprintedDbFailedMols.smi") #keeps track of failed molecules
    
    histogram = []
    
    for molecule in ifs.GetOEMols():
        mol, status = smiles2oemol(oechem.OEMolToSmiles(molecule))
        if (status):
            length = len([atom for atom in mol.GetAtoms()])
            histogram.append(length)
        else:
            oechem.OEWriteMolecule(ofs, mol)
    plt.hist(histogram, density=False, bins=40)
    plt.ylabel("# Of Molecules", fontsize=16)
    plt.xlabel("Atom Count", fontsize=16)
    plt.title("Fingerprinted Molecules: N-N Bonds", fontsize=24)
    plt.show()

def smiles2oemol(smiles):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)

    omega = oeomega.OEOmega()

    omega.SetMaxConfs(1)
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(True)
    omega.SetSampleHydrogens(True)
    omega.SetStrictStereo(True)
    omega.SetStrictAtomTypes(True)
    omega.SetIncludeInput(False)
    status = omega(mol)
    return (mol, status)

if __name__ == "__main__":
    makeHistogram()