"""
Various functions to create fingerprints for single nitrogen-nitrogen bonds
"""

from openeye import oechem, oeomega, oequacpac

# Object for calculating AM1 charges
AM1_CALCULATOR = oequacpac.OEAM1()


def find_neighboring_atoms(mol: oechem.OEMol, match: oechem.OEMatchBaseIter):
    """
    Finds and returns the neighbor atoms of a nitrogen-nitrogen bond as a list
    given a match. Number of neighboring atoms in the returned list can differ
    depending on how many neighbors the bond actually has and how bonds are
    found to be matching the SMIRKS pattern.
    """
    neighbors = []

    for trgtBond in match.Target().GetTargetBonds():
        trgtIdxs = (trgtBond.GetIdx() - 1, trgtBond.GetIdx())
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
    given a match. Number of Wiberg bond order values in the returned list can 
    differ depending on how many bonds are found to be matching the SMIRKS
    pattern.
    """
    wbo = []
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
    omega(mol)

    #Calculates bond orders of nitrogen-nitrogen bonds
    for trgtBond in match.Target().GetTargetBonds():
        trgtIdxs = (trgtBond.GetIdx() - 1, trgtBond.GetIdx())
        AM1_CALCULATOR.CalcAM1(results, mol)
        wbo.append(results.GetBondOrder(trgtIdxs[0], trgtIdxs[1]))

    return wbo
