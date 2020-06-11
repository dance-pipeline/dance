"""General testing utilities."""
import pathlib
from collections import Counter
from typing import List

from openeye import oechem

# pylint: disable=missing-function-docstring,invalid-name,unused-variable

#
# Functions
#


def oemol_from_smiles(smiles: str) -> oechem.OEMol:
    """oechem.OESmilesToMol wrapper that does not require creating a new OEMol."""
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    return mol


def get_canonical_isomeric_smiles(smiles: str) -> str:
    """Retrieve the canonical isomeric version of a SMILES string."""
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    return oechem.OEMolToSmiles(mol)


def get_list_of_canonical_isomeric_smiles(smiles_list: List[str]) -> List[str]:
    return [get_canonical_isomeric_smiles(smiles) for smiles in smiles_list]


def assert_smiles_in_file_are_equal(smiles_file: pathlib.Path, smiles: List[str]):
    """Checks that the SMILES in the given file are equivalent to those in the list (unordered)."""
    smiles_stream = oechem.oemolistream(str(smiles_file))
    outputted_smiles = Counter(oechem.OEMolToSmiles(mol) \
                            for mol in smiles_stream.GetOEMols())
    assert Counter(smiles) == outputted_smiles


def assert_ordered_smiles_in_file_are_equal(smiles_file: pathlib.Path, smiles: List[str]):
    """Checks that the SMILES in the given file are equivalent to those in the
    list."""
    smiles_stream = oechem.oemolistream(str(smiles_file))
    outputted_smiles = [oechem.OEMolToSmiles(mol) \
                            for mol in smiles_stream.GetOEMols()]
    assert smiles == outputted_smiles


def get_mols_from_oeb(oeb: pathlib.Path) -> List[oechem.OEMol]:
    """Returns an iterator over the molecules in the given OEB file."""
    assert oeb.is_file()  # Make sure the file exists.
    ifs = oechem.oemolistream(str(oeb))
    return ifs.GetOEMols()


def assert_smiles_in_oeb_are_equal(oeb: pathlib.Path, smiles: List[str]):
    """Checks that the molecules in the OEB have the same SMILES as those
    in the list."""
    outputted_smiles = Counter(oechem.OEMolToSmiles(mol) \
                            for mol in get_mols_from_oeb(oeb))
    assert Counter(smiles) == outputted_smiles


def assert_ordered_smiles_in_oeb_are_equal(oeb: pathlib.Path, smiles: List[str]):
    """Checks that the molecules in the OEB have the same SMILES as those in the
    list, and in the correct order."""
    outputted_smiles = [oechem.OEMolToSmiles(mol) for mol in get_mols_from_oeb(oeb)]
    assert smiles == outputted_smiles
