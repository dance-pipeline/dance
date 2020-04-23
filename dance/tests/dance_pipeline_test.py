"""Tests for DancePipeline."""
import pathlib
from typing import List

import pytest
from openeye import oechem

from dance import DancePipeline

# pylint: disable=missing-function-docstring,invalid-name

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


@pytest.fixture
def _pipeline_test_files(tmp_path):
    """Fixture for easily creating files for testing.

    Crates a SMILES file with the molecules from TEST_SMILES, an OEB filepath
    for the filter step, and an OEB filepath for the fingerprint step.
    """
    smiles_file = tmp_path / "smiles.smi"
    smiles_file.write_text("\n".join(TEST_SMILES))
    filter_output_oeb = tmp_path / "filter_output.oeb"
    fingerprint_output_oeb = tmp_path / "fingerprint_output.oeb"
    return smiles_file, filter_output_oeb, fingerprint_output_oeb


#
# Relevance functions
#

_NITROGEN = 7


def _relevant_always(mol: oechem.OEMol) -> bool:
    """Lets all molecules be marked as relevant."""
    # pylint: disable=unused-argument
    return True


def _relevant_if_contains_nitrogen(mol: oechem.OEMol) -> bool:
    return any(atom.GetAtomicNum() == _NITROGEN for atom in mol.GetAtoms())


#
# Test data
#

TEST_SMILES = ["N", "N#N", "O=C=O", "C#N"]
TEST_NUM_ATOMS_WITHOUT_EXPLICIT_HYDROGEN = [1, 2, 3, 2]  # Number of atoms in the molecules/SMILES above
TEST_OEMOLS = [_oemol_from_smiles(smiles) for smiles in TEST_SMILES]
TEST_CANONICAL_ISOMERIC_SMILES = _get_list_of_canonical_isomeric_smiles(TEST_SMILES)

#
# Initialization tests
#


def test_init_raises_exception_with_bad_database_type():
    with pytest.raises(RuntimeError):
        dp = DancePipeline("FOOBAR", "foobar.xyz")  # pylint: disable=unused-variable


def test_initial_molecules_not_available():
    dp = DancePipeline("SMILES", "foo.smi")
    assert dp.num_molecules is None


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
    # originally inputted, and the `num_molecules` attribute should have been
    # set correctly.
    assert dp.num_molecules == len(TEST_SMILES)
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
    assert dp.num_molecules == len(TEST_SMILES)
    _assert_smiles_in_oeb_are_equal(output_oeb, TEST_CANONICAL_ISOMERIC_SMILES)


def test_filters_molecules_with_relevance_function(tmp_path):
    smiles_file = tmp_path / "smiles.smi"
    smiles_file.write_text("\n".join(TEST_SMILES))
    output_oeb = tmp_path / "filter_output.oeb"

    for mol in TEST_OEMOLS:
        print(mol.NumAtoms())

    dp = DancePipeline("SMILES", smiles_file)
    dp.filter(_relevant_if_contains_nitrogen, output_oeb)

    assert dp.num_molecules == len(["N", "N#N", "C#N"])
    _assert_smiles_in_oeb_are_equal(output_oeb, \
        _get_list_of_canonical_isomeric_smiles(["N", "N#N", "C#N"]))


#
# Assign Fingerprint tests
#


@pytest.fixture
def _pipeline_executed_until_filter(_pipeline_test_files):
    """Provides a pipeline that has been executed up to and including the filter step.

    Also provides several associated files.
    """
    smiles_file, filter_output_oeb, fingerprint_output_oeb = _pipeline_test_files
    dp = DancePipeline("SMILES", smiles_file)
    dp.filter(_relevant_always, filter_output_oeb)
    return dp, smiles_file, filter_output_oeb, fingerprint_output_oeb


def test_assign_fingerprint_sets_output_oeb_attribute(_pipeline_executed_until_filter):
    dp, smiles_file, filter_output_oeb, fingerprint_output_oeb = _pipeline_executed_until_filter
    dp.assign_fingerprint(lambda mol: (), fingerprint_output_oeb)
    assert dp.fingerprint_output_oeb == fingerprint_output_oeb


def test_assigns_fingerprints_with_no_content(_pipeline_executed_until_filter):
    dp, smiles_file, filter_output_oeb, fingerprint_output_oeb = _pipeline_executed_until_filter

    dp.assign_fingerprint(lambda mol: (), fingerprint_output_oeb)

    # Check that the pipeline kept all the molecules after the fingerprint step.
    _assert_smiles_in_oeb_are_equal(fingerprint_output_oeb, TEST_CANONICAL_ISOMERIC_SMILES)

    # Check that the data tags on the molecules are correct.
    for mol in _get_mols_from_oeb(fingerprint_output_oeb):
        assert mol.GetIntData(dp.FINGERPRINT_LENGTH_NAME) == 0


def test_assigns_fingerprints_with_num_atoms_in_molecule(_pipeline_executed_until_filter):
    dp, smiles_file, filter_output_oeb, fingerprint_output_oeb = _pipeline_executed_until_filter

    dp.assign_fingerprint(lambda mol: (mol.NumAtoms(), ), fingerprint_output_oeb)

    _assert_smiles_in_oeb_are_equal(fingerprint_output_oeb, TEST_CANONICAL_ISOMERIC_SMILES)
    for mol, num_atoms in zip(_get_mols_from_oeb(fingerprint_output_oeb), TEST_NUM_ATOMS_WITHOUT_EXPLICIT_HYDROGEN):
        assert mol.GetIntData(dp.FINGERPRINT_LENGTH_NAME) == 1
        assert mol.GetDoubleData(f"{dp.FINGERPRINT_VALUE_NAME}_0") == num_atoms


def test_assigns_fingerprints_with_multiple_entries(_pipeline_executed_until_filter):
    dp, smiles_file, filter_output_oeb, fingerprint_output_oeb = _pipeline_executed_until_filter

    # The fingerprint is based on the number of atoms, along with a value of 100
    # at the end.
    dp.assign_fingerprint(lambda mol: (mol.NumAtoms(), mol.NumAtoms() - 1, mol.NumAtoms() - 2, 100),
                          fingerprint_output_oeb)

    _assert_smiles_in_oeb_are_equal(fingerprint_output_oeb, TEST_CANONICAL_ISOMERIC_SMILES)
    for mol, num_atoms in zip(_get_mols_from_oeb(fingerprint_output_oeb), TEST_NUM_ATOMS_WITHOUT_EXPLICIT_HYDROGEN):
        assert mol.GetIntData(dp.FINGERPRINT_LENGTH_NAME) == 4
        for i in range(3):
            assert mol.GetDoubleData(f"{dp.FINGERPRINT_VALUE_NAME}_{i}") == num_atoms - i
        assert mol.GetDoubleData(f"{dp.FINGERPRINT_VALUE_NAME}_3") == 100
