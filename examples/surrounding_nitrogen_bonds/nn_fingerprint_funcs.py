"""
Various functions to create fingerprints for single nitrogen-nitrogen bonds
"""

from openeye import oechem, oeomega, oequacpac
import logging

# Configure logger.
logger = logging.getLogger(__name__)

# Object for calculating AM1 charges
AM1_CALCULATOR = oequacpac.OEAM1()

def smiles2oemol(smiles):
    """
    Conforms molecules from a smiles string to an OEMol with all the traits
    necessary for OpenEye functions to return accurate data
    """
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    #initialize omega
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(True)
    omega.SetSampleHydrogens(True)
    omega.SetStrictStereo(False)
    omega.SetStrictAtomTypes(True)
    omega.SetIncludeInput(False)
    status = omega(mol)
    
    return (mol, status)

def find_neighboring_atoms(mol: oechem.OEMol, match: oechem.OEMatchBaseIter):
    """
    Finds and returns the atomic numbers of neighboring atoms to the 
    nitrogen-nitrogen bond of a molecule as a list. The maximum number of 
    potential neighbors to be used is set to 4. Molecules that do not have
    4 neighbors will have a 0 in the remaining spots. 
    """
    neighbors = [0, 0, 0, 0]
    indx = 0
    
    target_bond = match.Target().GetTargetBonds().Target()
    target_atoms = (target_bond.GetBgn(), target_bond.GetEnd())
    
    for bond in mol.GetBonds():
        if(bond.GetBgn() in target_atoms or bond.GetEnd() in target_atoms):
            if bond.GetBgn() not in target_atoms:
                neighbors[indx] = bond.GetBgn().GetAtomicNum()
                indx += 1
            elif bond.GetEnd() not in target_atoms:
                neighbors[indx] = bond.GetEnd().GetAtomicNum()
                indx += 1
            if indx == 4:
                break
    
    return neighbors

def wiberg_bond_order(mol: oechem.OEMol, match: oechem.OEMatchBaseIter):
    """
    Finds and returns the Wiberg bond order of the nitrogen-nitrogen bond 
    of a molecule.
    """
    results = oequacpac.OEAM1Results()

    target_bond = match.Target().GetTargetBonds().Target()
    target_idxs = (target_bond.GetBgn().GetIdx(), target_bond.GetEnd().GetIdx())
    
    #Calculates bond orders of nitrogen-nitrogen bonds
    AM1_CALCULATOR.CalcAM1(results, mol)
    
    return results.GetBondOrder(target_idxs[0], target_idxs[1])
