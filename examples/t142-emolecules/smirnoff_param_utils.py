"""Utilities for scripts working with SMIRNOFF params."""
import json
import logging
from collections import defaultdict
from typing import Dict, List

from openeye import oechem, oeomega, oequacpac
from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField

# Configure logger.
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name

# smirnoff99Frosst force field
FORCE_FIELD = ForceField('test_forcefields/smirnoff99Frosst.offxml')

# Name of the tag that holds the parameters in the molecule. Parameters
# are simply stored as JSON strings.
PARAM_TAG_NAME = "SMIRNOFF_PARAMS"

# Object for calculating AM1 charges - instead of reconstructing everywhere
AM1_CALCULATOR = oequacpac.OEAM1()


def write_params_to_mol(mol: oechem.OEMol, params: Dict[str, List[List[int]]]):
    """Writes the given list of parameters to tags in the given molecule."""
    mol.SetStringData(PARAM_TAG_NAME, json.dumps(params))


def read_params_from_mol(mol: oechem.OEMol) -> Dict[str, List[List[int]]]:
    """Reads the parameters from the tags in the given molecule."""
    return json.loads(mol.GetStringData(PARAM_TAG_NAME))


def calculate_mol_params(mol: oechem.OEMol) -> Dict[str, List[List[int]]]:
    """Calculates parameters of the given molecule.

    Returns a dict where the keys are parameter ids and the values are lists
    of indices where the parameter occurs (each entry in the list is itself a
    list because the parameter involves multiple atoms).
    """
    oechem.OEAddExplicitHydrogens(mol)
    off_mol = Molecule.from_openeye(mol, allow_undefined_stereo=True)
    topology = Topology.from_molecules(off_mol)
    molecule_force_list = FORCE_FIELD.label_molecules(topology)

    params = defaultdict(list)
    for _, force_dict in molecule_force_list[0].items():
        for (atom_indices, parameter) in force_dict.items():
            params[parameter.id].append(atom_indices)

    return params


def count_total_molecules(filename: str) -> int:
    """Counts number of molecules in a file"""
    total = 0
    ifs = oechem.oemolistream(filename)
    mol = oechem.OEMol()
    while oechem.OEReadMolecule(ifs, mol):
        total += 1
    return total


def calculate_t142_central_wbo(mol: oechem.OEMol, params: Dict[str, List[List[int]]]) -> float:
    """Calculates the WBO between the central atoms in the t142 param in the molecule.

    (WBO is Wiberg Bond Order.)

    The `params` argument contains the parameters of the molecule (see
    `calculate_mol_params`).

    Returns -1 if the calculation fails.
    """
    # Only use first occurrence of the parameter.
    indices = params['t142'][0]

    # For torsion parameters such as t142, the central atoms should be at the
    # second and third index.
    central_indices = [indices[1], indices[2]]

    # Generate molecule conformer.
    oechem.OEAddExplicitHydrogens(mol)
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetCanonOrder(False)
    omega.SetSampleHydrogens(True)
    omega.SetEnergyWindow(15.0)  #unit?
    omega.SetRMSThreshold(1.0)
    # Don't generate random stereoisomer if not specified.
    omega.SetStrictStereo(True)
    status = omega(mol)

    if status is False:
        omega.SetStrictStereo(False)
        new_status = omega(mol)
        if new_status is False:
            logger.error("Failed to generate conformer for %s", oechem.OEMolToSmiles(mol))
            return -1

    # Calculate the WBO between the two central atoms.
    conf = next(iter(mol.GetConfs()))
    charged_copy = oechem.OEMol(conf)
    results = oequacpac.OEAM1Results()

    if not AM1_CALCULATOR.CalcAM1(results, charged_copy):
        logger.error("Failed to assign partial charges to %s", oechem.OEMolToSmiles(mol))
        return -2

    return results.GetBondOrder(central_indices[0], central_indices[1])
