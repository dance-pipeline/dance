"""
Script for generating a histogram based on the bond length of
single nitrogen-nitrogen bonds
Usage:
    python nn_bond_length_histogram.py --smiles_database [DATABASE] --fingerprinted [T or F]
"""

from openeye import oechem, oeomega, oequacpac, oeff
import argparse
import matplotlib.pyplot as plt
import math

def main():
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
    
    generate_histogram(args)

def generate_histogram(args):
    """
    Generates a bond length histogram for single nitrogen-nitrogen bonds
    using a .smi file
    """
    ifs = oechem.oemolistream(args.smiles_database)
    subs = oechem.OESubSearch("[#7:1]-[#7:2]")
    bond_lengths = []
    
    for molecule in ifs.GetOEMols():
        mol, status = smiles2oemol(oechem.OEMolToSmiles(molecule))
        if status:
            optimize_mol(mol)
            
            match = subs.Match(mol, True)
            target_bond = match.Target().GetTargetBonds().Target()
            target_idxs = (target_bond.GetBgn().GetIdx(), target_bond.GetEnd().GetIdx())
            
            n1_coords = mol.GetCoords()[target_idxs[0]]
            n2_coords = mol.GetCoords()[target_idxs[1]]
            
            bond_length = get_bond_length(n1_coords, n2_coords)
            bond_lengths.append(bond_length)
    
    plt.hist(bond_lengths, density=True, bins=30)
    
    if args.fingerprinted.lower() == "t":
        plt.title(f"Fingerprinted N-N Molecules ({len(bond_lengths)} total mols)")
    elif args.fingerprinted.lower() == "f":
        plt.title(f"Random N-N Molecules ({len(bond_lengths)} total mols)")
    plt.ylabel("Frequency")
    plt.xlabel("Bond Length (angstrom)")
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
    
def optimize_mol(mol):
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

def get_bond_length(n1_coords, n2_coords):
    """
    Calculates and returns the bond length between 2 nitrogen atoms
    """
    return math.sqrt( (n2_coords[0] - n1_coords[0])**2 +
                     (n2_coords[1] - n1_coords[1])**2 +
                     (n2_coords[2] - n1_coords[2])**2 )

if __name__ == "__main__":
    main()
    
