"""Tests for DancePipeline.
"""
from dance import DancePipeline
import pathlib
import pytest
from openeye import oechem
from typing import List

#
# Utilities
#


def _oemol_from_smiles(smiles: str) -> oechem.OEMol:
    """oechem.OESmilesToMol wrapper that does not require creating a new OEMol."""
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    return mol


def _get_canonical_isomeric_smiles(smiles: str) -> str:
    """Retrieve the canonical isomeric version of a SMILES string."""
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    return oechem.OEMolToSmiles(mol)


def _get_list_of_canonical_isomeric_smiles(smiles_list: List[str]) -> List[str]:
    return [_get_canonical_isomeric_smiles(smiles) for smiles in smiles_list]


def _get_mols_from_oeb(oeb: pathlib.Path) -> List[oechem.OEMol]:
    """Returns an iterator over the molecules in the given OEB file."""
    assert oeb.is_file()  # Make sure the file exists.
    ifs = oechem.oemolistream(str(oeb))
    return ifs.GetOEMols()


def _assert_smiles_in_oeb_are_equal(oeb: pathlib.Path, smiles: List[str]):
    """Checks that the molecules in the oeb have the same SMILES as those in the list."""
    outputted_smiles = set(oechem.OEMolToSmiles(mol) \
                            for mol in _get_mols_from_oeb(oeb))
    assert set(smiles) == outputted_smiles


#
# Relevance functions
#

_NITROGEN = 7


def _relevant_always(mol: oechem.OEMol) -> bool:
    """Lets all molecules be marked as relevant."""
    return True


def _relevant_if_contains_nitrogen(mol: oechem.OEMol) -> bool:
    return any(atom.GetAtomicNum() == _NITROGEN for atom in mol.GetAtoms())


#
# Test data
#

TEST_SMILES = ["N", "N#N", "O=C=O", "C#N"]
TEST_OEMOLS = [_oemol_from_smiles(smiles) for smiles in TEST_SMILES]
TEST_CANONICAL_ISOMERIC_SMILES = _get_list_of_canonical_isomeric_smiles(TEST_SMILES)

#
# Initialization tests
#


def test_init_raises_exception_with_bad_database_type():
    with pytest.raises(RuntimeError):
        dp = DancePipeline("FOOBAR", "foobar.xyz")


#
# Filter tests
#


def test_filter_sets_output_oeb_attribute(tmp_path):
    smiles_file = tmp_path / "smiles.smi"
    smiles_file.write_text("\n".join(TEST_SMILES))
    output_oeb = tmp_path / "filter_output.oeb"

    dp = DancePipeline("SMILES", smiles_file)
    dp.filter(_relevant_always, output_oeb)

    assert dp.filter_output_oeb == output_oeb


def test_filters_from_smiles_database(tmp_path):
    smiles_file = tmp_path / "smiles.smi"
    smiles_file.write_text("\n".join(TEST_SMILES))
    output_oeb = tmp_path / "filter_output.oeb"

    dp = DancePipeline("SMILES", smiles_file)
    dp.filter(_relevant_always, output_oeb)

    # After "filtering," the output oeb should have all the molecules that were
    # originally inputted.
    _assert_smiles_in_oeb_are_equal(output_oeb, TEST_CANONICAL_ISOMERIC_SMILES)


def test_filters_from_mol2_dir_database(tmp_path):
    mol2dir = tmp_path / "molecules"
    mol2dir.mkdir()
    for idx, mol in enumerate(TEST_OEMOLS):
        mol2file = mol2dir / f"{idx}.mol2"
        ofs = oechem.oemolostream(str(mol2file))
        oechem.OEWriteMolecule(ofs, mol)
        ofs.close()
    output_oeb = tmp_path / "filter_output.oeb"

    dp = DancePipeline("MOL2_DIR", mol2dir)
    dp.filter(_relevant_always, output_oeb)

    # Even though this is a mol2dir, we still use SMILES to make sure the
    # molecules are equal, as we do not have a way to compare two OEMols.
    _assert_smiles_in_oeb_are_equal(output_oeb, TEST_CANONICAL_ISOMERIC_SMILES)


def test_filters_molecules_with_relevance_function(tmp_path):
    smiles_file = tmp_path / "smiles.smi"
    smiles_file.write_text("\n".join(TEST_SMILES))
    output_oeb = tmp_path / "filter_output.oeb"

    for mol in TEST_OEMOLS:
        print(mol.NumAtoms())

    dp = DancePipeline("SMILES", smiles_file)
    dp.filter(_relevant_if_contains_nitrogen, output_oeb)

    _assert_smiles_in_oeb_are_equal(output_oeb, \
        _get_list_of_canonical_isomeric_smiles(["N", "N#N", "C#N"]))
