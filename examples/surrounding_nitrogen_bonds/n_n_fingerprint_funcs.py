"""
Various functions to create fingerprints for single nitrogen-nitrogen bonds
"""

from openeye import oechem, oeomega, oequacpac

def find_neighboring_atoms(mol: oechem.OEMol, match: oechem.OEMatchBaseIter):
    """
    Finds and returns the neighbor atoms of a nitrogen-nitrogen bond as a list
    given a match
    """
    neighbors = []
    
    for trgtBond in match.Target().GetTargetBonds():
        trgtIdxs = (trgtBond.GetIdx()-1, trgtBond.GetIdx())
        for atom in mol.GetAtoms():
            if atom.GetIdx() in trgtIdxs:
                for bond in atom.GetBonds():
                    nbor = bond.GetNbr(atom)
                    if nbor.GetIdx() not in trgtIdxs:
                        neighbors.append(nbor.GetAtomicNum())
    return neighbors

def wiberg_bond_order(mol: oechem.OEMol, match: oechem.OEMatchBaseIter):
    """
    Finds and returns a list of Wiberg bond orders of nitrogen-nitrogen bonds
    given a match
    """
    wbo = []
    AM1_CALCULATOR = oequacpac.OEAM1()
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
    
    #Calculates bond orders of nitrogen-nitrogen bonds
    for trgtBond in match.Target().GetTargetBonds():
        trgtIdxs = (trgtBond.GetIdx()-1, trgtBond.GetIdx())
        AM1_CALCULATOR.CalcAM1(results, mol)
        wbo.append(results.GetBondOrder(trgtIdxs[0], trgtIdxs[1]))
        
    return wbo