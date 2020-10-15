"""
Various functions to create fingerprints for single nitrogen-nitrogen bonds
"""

from openeye import oechem, oeomega, oequacpac
import logging

# Configure logger.
logger = logging.getLogger(__name__)

# Object for calculating AM1 charges
AM1_CALCULATOR = oequacpac.OEAM1()


def find_neighboring_atoms(mol: oechem.OEMol, match: oechem.OEMatchBaseIter):
    """
    Finds and returns the atomic numbers of neighboring atoms to the 
    nitrogen-nitrogen bond of a molecule as a list. The maximum number of 
    potential neighbors to be used is set to 3.
    """
    neighbors = [0 for _ in range(3)]
    indx = 0

    trgtBond = match.Target().GetTargetBonds().Target()
    trgtIdxs = (trgtBond.GetIdx() - 1, trgtBond.GetIdx())
    for atom in mol.GetAtoms():
        if atom.GetIdx() in trgtIdxs:
            for bond in atom.GetBonds():
                nbor = bond.GetNbr(atom)
                if nbor.GetIdx() not in trgtIdxs:
                    neighbors[indx] = nbor.GetAtomicNum()
                    indx += 1
                    if indx == 3: 
                        return neighbors
    return neighbors


def wiberg_bond_order(mol: oechem.OEMol, match: oechem.OEMatchBaseIter):
    """
    Finds and returns the Wiberg bond order of the nitrogen-nitrogen bond 
    of a molecule.
    Returns -1 if the calculation fails.
    """
    results = oequacpac.OEAM1Results()

    #Generates molecule conformer
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(True)
    omega.SetSampleHydrogens(True)
    omega.SetStrictStereo(True)
    omega.SetStrictAtomTypes(True)
    omega.SetIncludeInput(False)
    status = omega(mol)
    
    if status is False:
        omega.SetStrictStereo(False)
        newStatus = omega(mol)
        if newStatus is False:
            logger.error("Failed to generate conformer for %s", oechem.OEMolToSmiles(mol))
            return [-1]

    #Calculates bond orders of nitrogen-nitrogen bonds
    trgtBond = match.Target().GetTargetBonds().Target()
    trgtIdxs = (trgtBond.GetIdx() - 1, trgtBond.GetIdx())
    AM1_CALCULATOR.CalcAM1(results, mol)
    return results.GetBondOrder(trgtIdxs[0], trgtIdxs[1])
    
