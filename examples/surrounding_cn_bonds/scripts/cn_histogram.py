"""
Script for generating a histogram of the C-N molecule substructures
in random/fingerprinted sets of molecules. 
Usage:
    python cn_histogram.py --smiles_database [DATABASE] --fingerprinted [T or F]
"""

from openeye import oechem, oeomega, oequacpac, oeff
import argparse
import matplotlib.pyplot as plt
import math

def generateHistogram():
    """
    Generates a bond length histogram for single nitrogen-nitrogen bonds
    using a .smi file
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles_database",
                        type=str,
                        required=True,
                        help=("Name of a SMILES file containing molecules "
                              "to generate a bond length histogram."))
    parser.add_argument("--fingerprinted",
                        type=str,
                        required=True,
                        help=("True or False if the database is fingerprinted or not"))
    args = parser.parse_args()
    
    ifs = oechem.oemolistream(args.smiles_database)
    subs = oechem.OESubSearch("[#6:1]-[#7:2]")
    histogram = []
    mol_count = 0
    for molecule in ifs.GetOEMols():
        mol, status = smiles2oemol(oechem.OEMolToSmiles(molecule))
        if (status):
            if subs.SingleMatch(mol):
                mol_count = mol_count + 1
                length = len([atom for atom in mol.GetAtoms()])
                histogram.append(length)
    title = ""
    colour = "blue"
    if (args.fingerprinted.lower() == "t"):
        title = f"Atom Count of {mol_count} Fingerprinted C-N Molecules"
    elif (args.fingerprinted.lower() == "f"):
        title = f"Atom Count of {mol_count} Random C-N Molecules"
        colour = "red"
    plt.hist(histogram, range=[0,100],color = colour, density=False, bins=30)
    plt.ylim(0,100)
    plt.title(title)
    plt.ylabel("Frequency")
    plt.xlabel("Atom Count")
    plt.show(block=True)

def smiles2oemol(smiles):
    """
    Conforms molecules from a smiles string to an OEMol with all the traits
    necessary for OpenEye functions to return accurate data
    """
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    # Initialize Omega
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(True)
    omega.SetSampleHydrogens(True)  # Word to the wise: skipping this step can lead to significantly different charges!
    omega.SetStrictStereo(True)
    omega.SetStrictAtomTypes(True)
    omega.SetIncludeInput(False) # don't include input
    status = omega(mol)
    return (mol, status)
    
def optimizeMol(mol):
    """
    Optimizes OEMols using an AM1 calculation so that they return the most 
    accurate data possible 
    """
    am1 = oequacpac.OEAM1()  # this is a low-level api class that will not perform any optimization on its own
    results = oequacpac.OEAM1Results()  # this is where you can access the energy as well as the bond orders
    vecCoords = oechem.OEDoubleArray(3 * mol.NumAtoms())
    
    for conf in mol.GetConfs():
        conf.GetCoords(vecCoords)
        am1.CalcAM1(results, conf)
        optimizer = oeff.OEBFGSOpt()
        optimizer(am1, vecCoords, vecCoords)  # (MolFunc1, Input Coords, Output Coords)
        conf.SetCoords(vecCoords)
        am1.CalcAM1(results, conf)

    return mol   

if __name__ == "__main__":
    generateHistogram()