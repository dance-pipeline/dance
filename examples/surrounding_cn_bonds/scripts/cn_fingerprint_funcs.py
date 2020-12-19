from openeye import oechem, oeomega, oequacpac
import logging

# Configure logger.
logger = logging.getLogger(__name__)

# Object for calculating AM1 charges
AM1_CALCULATOR = oequacpac.OEAM1()

def smiles_to_oemol(smiles):
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
    
    if status is False:
        omega.SetStrictStereo(False)
        newStatus = omega(mol)
        if newStatus is False:
            logger.error("Failed to generate conformer for %s", oechem.OEMolToSmiles(mol))
            
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
    
    targetBond = match.Target().GetTargetBonds().Target()
    targetAtoms = (targetBond.GetBgn(), targetBond.GetEnd())
    
    for bond in mol.GetBonds():
        if(bond.GetBgn() in targetAtoms or bond.GetEnd() in targetAtoms):
            if bond.GetBgn() not in targetAtoms:
                neighbors[indx] = bond.GetBgn().GetAtomicNum()
                indx += 1
            elif bond.GetEnd() not in targetAtoms:
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

    targetBond = match.Target().GetTargetBonds().Target()
    targetIdxs = (targetBond.GetBgn().GetIdx(), targetBond.GetEnd().GetIdx())
    
    #Calculates bond orders of nitrogen-nitrogen bonds
    AM1_CALCULATOR.CalcAM1(results, mol)
    
    return results.GetBondOrder(targetIdxs[0], targetIdxs[1])