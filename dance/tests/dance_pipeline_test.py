"""Tests for DancePipeline."""
import random

import pytest
from openeye import oechem

from dance import DancePipeline

from . import utils

# pylint: disable=missing-function-docstring,invalid-name,unused-variable

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
TEST_OEMOLS = [utils.oemol_from_smiles(smiles) for smiles in TEST_SMILES]
TEST_CANONICAL_ISOMERIC_SMILES = utils.get_list_of_canonical_isomeric_smiles(TEST_SMILES)

#
# Initialization tests
#


def test_init_raises_exception_with_bad_database_type():
    with pytest.raises(RuntimeError):
        dp = DancePipeline("FOOBAR", "foobar.baz")


def test_initial_molecules_not_available():
    dp = DancePipeline("SMILES", "foo.smi")
    assert dp.num_molecules is None


#
# Utility tests
#


def test_fingerprint_retrieval_raises_exceptions_with_bad_data():
    # Retrieving from a molecule with no length tag.
    mol = oechem.OEMol()
    with pytest.raises(ValueError):
        DancePipeline.get_fingerprint_from_mol(mol)

    # Retrieving from a molecule with length tag but missing values.
    mol.SetIntData(DancePipeline.FINGERPRINT_LENGTH_NAME, 4)
    with pytest.raises(ValueError):
        DancePipeline.get_fingerprint_from_mol(mol)


def test_retrieves_correct_fingerprint_from_mol():
    mol = oechem.OEMol()
    mol.SetIntData(DancePipeline.FINGERPRINT_LENGTH_NAME, 4)
    for i, val in enumerate([2, 1, 8, 7]):
        mol.SetDoubleData(f"{DancePipeline.FINGERPRINT_VALUE_NAME}_{i}", val)

    assert DancePipeline.get_fingerprint_from_mol(mol) == (2, 1, 8, 7)


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
    utils.assert_smiles_in_oeb_are_equal(output_oeb, TEST_CANONICAL_ISOMERIC_SMILES)


def test_filters_from_sdf_database(tmp_path):
    # Create the SDF file and write the molecules to it.
    sdf_file = tmp_path / "molecules.sdf"
    sdf_stream = oechem.oemolostream(str(sdf_file))
    for smiles in TEST_SMILES:
        mol = utils.oemol_from_smiles(smiles)
        oechem.OEWriteMolecule(sdf_stream, mol)
    sdf_stream.close()

    output_oeb = tmp_path / "filter_output.oeb"

    dp = DancePipeline("SDF", sdf_file)
    dp.filter(_relevant_always, output_oeb)

    assert dp.num_molecules == len(TEST_SMILES)
    utils.assert_smiles_in_oeb_are_equal(output_oeb, TEST_CANONICAL_ISOMERIC_SMILES)


def test_filters_from_oeb_database(tmp_path):
    # Create the OEB file and write the molecules to it.
    oeb_file = tmp_path / "molecules.oeb"
    oeb_stream = oechem.oemolostream(str(oeb_file))
    for smiles in TEST_SMILES:
        mol = utils.oemol_from_smiles(smiles)
        oechem.OEWriteMolecule(oeb_stream, mol)
    oeb_stream.close()

    output_oeb = tmp_path / "filter_output.oeb"

    dp = DancePipeline("OEB", oeb_file)
    dp.filter(_relevant_always, output_oeb)

    assert dp.num_molecules == len(TEST_SMILES)
    utils.assert_smiles_in_oeb_are_equal(output_oeb, TEST_CANONICAL_ISOMERIC_SMILES)


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
    utils.assert_smiles_in_oeb_are_equal(output_oeb, TEST_CANONICAL_ISOMERIC_SMILES)


def test_filters_molecules_with_relevance_function(tmp_path):
    smiles_file = tmp_path / "smiles.smi"
    smiles_file.write_text("\n".join(TEST_SMILES))
    output_oeb = tmp_path / "filter_output.oeb"

    dp = DancePipeline("SMILES", smiles_file)
    dp.filter(_relevant_if_contains_nitrogen, output_oeb)

    assert dp.num_molecules == len(["N", "N#N", "C#N"])
    utils.assert_smiles_in_oeb_are_equal(output_oeb, \
        utils.get_list_of_canonical_isomeric_smiles(["N", "N#N", "C#N"]))


#
# Assign Fingerprint tests
#


@pytest.fixture
def _pipeline_executed_until_filter(pipeline_test_files):
    """Provides a pipeline that has been executed up to and including the filter step.

    Also provides several associated files.
    """
    (smiles_file, filter_output_oeb, fingerprint_output_oeb, smiles_dataset_file,
     sorted_by_fingerprint_oeb) = pipeline_test_files

    smiles_file.write_text("\n".join(TEST_SMILES))

    dp = DancePipeline("SMILES", smiles_file)
    dp.filter(_relevant_always, filter_output_oeb)

    return dp, smiles_file, filter_output_oeb, fingerprint_output_oeb


def test_assign_fingerprint_sets_output_oeb_attribute(_pipeline_executed_until_filter):
    (dp, smiles_file, filter_output_oeb, fingerprint_output_oeb) = _pipeline_executed_until_filter

    dp.assign_fingerprint(lambda mol: (), fingerprint_output_oeb)

    assert dp.fingerprint_output_oeb == fingerprint_output_oeb


def test_assigns_fingerprints_with_no_content(_pipeline_executed_until_filter):
    (dp, smiles_file, filter_output_oeb, fingerprint_output_oeb) = _pipeline_executed_until_filter

    dp.assign_fingerprint(lambda mol: (), fingerprint_output_oeb)

    # Check that the pipeline kept all the molecules after the fingerprint step.
    utils.assert_smiles_in_oeb_are_equal(fingerprint_output_oeb, TEST_CANONICAL_ISOMERIC_SMILES)

    # Check that the data tags on the molecules are correct.
    for mol in utils.get_mols_from_oeb(fingerprint_output_oeb):
        assert mol.GetIntData(dp.FINGERPRINT_LENGTH_NAME) == 0


def test_assigns_fingerprints_with_num_atoms_in_molecule(_pipeline_executed_until_filter):
    (dp, smiles_file, filter_output_oeb, fingerprint_output_oeb) = _pipeline_executed_until_filter

    dp.assign_fingerprint(lambda mol: (mol.NumAtoms(), ), fingerprint_output_oeb)

    utils.assert_smiles_in_oeb_are_equal(fingerprint_output_oeb, TEST_CANONICAL_ISOMERIC_SMILES)
    for mol, num_atoms in zip(utils.get_mols_from_oeb(fingerprint_output_oeb),
                              TEST_NUM_ATOMS_WITHOUT_EXPLICIT_HYDROGEN):
        assert mol.GetIntData(dp.FINGERPRINT_LENGTH_NAME) == 1
        assert mol.GetDoubleData(f"{dp.FINGERPRINT_VALUE_NAME}_0") == num_atoms


def test_assigns_fingerprints_with_multiple_entries(_pipeline_executed_until_filter):
    (dp, smiles_file, filter_output_oeb, fingerprint_output_oeb) = _pipeline_executed_until_filter

    # The fingerprint is based on the number of atoms, along with a value of 100
    # at the end.
    dp.assign_fingerprint(lambda mol: (mol.NumAtoms(), mol.NumAtoms() - 1, mol.NumAtoms() - 2, 100),
                          fingerprint_output_oeb)

    utils.assert_smiles_in_oeb_are_equal(fingerprint_output_oeb, TEST_CANONICAL_ISOMERIC_SMILES)
    for mol, num_atoms in zip(utils.get_mols_from_oeb(fingerprint_output_oeb),
                              TEST_NUM_ATOMS_WITHOUT_EXPLICIT_HYDROGEN):
        assert mol.GetIntData(dp.FINGERPRINT_LENGTH_NAME) == 4
        for i in range(3):
            assert mol.GetDoubleData(f"{dp.FINGERPRINT_VALUE_NAME}_{i}") == num_atoms - i
        assert mol.GetDoubleData(f"{dp.FINGERPRINT_VALUE_NAME}_3") == 100


#
# Selection tests
#


@pytest.fixture
def _pipeline_executed_until_fingerprint(pipeline_test_files):
    """Provides a pipeline that has been executed up to and including the assign_fingerprint step.

    The fingerprint function assigns a fingerprint consisting of the number of
    atoms in the molecule.

    Also provides several associated files.
    """
    (smiles_file, filter_output_oeb, fingerprint_output_oeb, smiles_dataset_file,
     sorted_by_fingerprint_oeb) = pipeline_test_files

    smiles_file.write_text("\n".join(TEST_SMILES))

    dp = DancePipeline("SMILES", smiles_file)
    dp.filter(_relevant_always, filter_output_oeb)
    dp.assign_fingerprint(lambda mol: (mol.NumAtoms(), ), fingerprint_output_oeb)

    return dp, smiles_file, filter_output_oeb, fingerprint_output_oeb, smiles_dataset_file, sorted_by_fingerprint_oeb


def test_select_raises_exception_with_bad_output_type(_pipeline_executed_until_fingerprint):
    (dp, *unused) = _pipeline_executed_until_fingerprint
    with pytest.raises(RuntimeError):
        dp.select(1, "FOOBAR", "foobar.baz")


def test_select_raises_exception_with_bad_selection_freq(_pipeline_executed_until_fingerprint):
    (dp, *unused, smiles_dataset_file, unused2) = _pipeline_executed_until_fingerprint

    with pytest.raises(RuntimeError):
        dp.select(0, "SMILES", smiles_dataset_file)


def test_select_sets_sorted_fingerprint_oeb_attribute(_pipeline_executed_until_fingerprint):
    (dp, *unused, smiles_dataset_file, sorted_by_fingerprint_oeb) = _pipeline_executed_until_fingerprint

    dp.select(1, "SMILES", smiles_dataset_file, sorted_by_fingerprint_oeb)

    assert dp.sorted_by_fingerprint_oeb == sorted_by_fingerprint_oeb


def test_select_outputs_to_oeb(_pipeline_executed_until_fingerprint, tmp_path):
    (dp, *unused, smiles_dataset_file, sorted_by_fingerprint_oeb) = _pipeline_executed_until_fingerprint
    oeb_dataset = tmp_path / "dataset.oeb"

    dp.select(1, "OEB", oeb_dataset, sorted_by_fingerprint_oeb)

    utils.assert_smiles_in_oeb_are_equal(oeb_dataset, TEST_CANONICAL_ISOMERIC_SMILES)


def test_select_can_choose_every_molecule(_pipeline_executed_until_fingerprint):
    (dp, *unused, smiles_dataset_file, sorted_by_fingerprint_oeb) = _pipeline_executed_until_fingerprint
    dp.select(1, "SMILES", smiles_dataset_file, sorted_by_fingerprint_oeb)

    # Either ordering is okay for the sorted output, as both C#N and N#N have
    # two atoms (not counting explicit hydrogens, that is).
    outputted_smiles = \
            [oechem.OEMolToSmiles(mol) for mol in utils.get_mols_from_oeb(sorted_by_fingerprint_oeb)]
    assert outputted_smiles == \
            utils.get_list_of_canonical_isomeric_smiles(["N", "N#N", "C#N", "O=C=O"]) or \
           outputted_smiles == \
            utils.get_list_of_canonical_isomeric_smiles(["N", "C#N", "N#N", "O=C=O"])

    utils.assert_smiles_in_file_are_equal(smiles_dataset_file, TEST_CANONICAL_ISOMERIC_SMILES)


def test_select_only_certain_molecules(_pipeline_executed_until_fingerprint):
    (dp, *unused, smiles_dataset_file, sorted_by_fingerprint_oeb) = _pipeline_executed_until_fingerprint

    # Setting the threshold to 3 ensures we use multiple files in the sorting.
    dp.select(3, "SMILES", smiles_dataset_file, sorted_by_fingerprint_oeb, in_memory_sorting_threshold=3)

    utils.assert_smiles_in_file_are_equal(smiles_dataset_file,
                                          utils.get_list_of_canonical_isomeric_smiles(["N", "O=C=O"]))


def test_select_with_longer_fingerprints(pipeline_test_files):
    (smiles_file, filter_output_oeb, fingerprint_output_oeb, smiles_dataset_file,
     sorted_by_fingerprint_oeb) = pipeline_test_files

    smiles_file.write_text("\n".join(TEST_SMILES))

    dp = DancePipeline("SMILES", smiles_file)
    dp.filter(_relevant_always, filter_output_oeb)

    # Each fingerprint consists of (1, 1, num_atoms), to ensure that fingerprints
    # are sorted correctly when the fingerprint is longer.
    dp.assign_fingerprint(lambda mol: (1, 1, mol.NumAtoms()), fingerprint_output_oeb)
    dp.select(3, "SMILES", smiles_dataset_file, sorted_by_fingerprint_oeb, in_memory_sorting_threshold=3)

    # Check that the molecules are sorted by fingerprint.
    outputted_smiles = \
            [oechem.OEMolToSmiles(mol) for mol in utils.get_mols_from_oeb(sorted_by_fingerprint_oeb)]
    assert outputted_smiles == \
            utils.get_list_of_canonical_isomeric_smiles(["N", "N#N", "C#N", "O=C=O"]) or \
           outputted_smiles == \
            utils.get_list_of_canonical_isomeric_smiles(["N", "C#N", "N#N", "O=C=O"])

    # Check that the correct molecules were selected.
    utils.assert_smiles_in_file_are_equal(smiles_dataset_file,
                                          utils.get_list_of_canonical_isomeric_smiles(["N", "O=C=O"]))


def test_select_on_larger_dataset(pipeline_test_files):
    (smiles_file, filter_output_oeb, fingerprint_output_oeb, smiles_dataset_file,
     sorted_by_fingerprint_oeb) = pipeline_test_files

    # The pipeline does not filter out repeated molecules, so this is okay.
    smiles = ["N", "O=C=O", "C#N"] * 10
    random.shuffle(smiles)  # The pipeline should work for any molecule ordering.
    smiles_file.write_text("\n".join(smiles))

    # Pipeline execution.
    dp = DancePipeline("SMILES", smiles_file)
    dp.filter(_relevant_always, filter_output_oeb)
    dp.assign_fingerprint(lambda mol: (mol.NumAtoms(), ), fingerprint_output_oeb)
    dp.select(10, "SMILES", smiles_dataset_file, sorted_by_fingerprint_oeb, in_memory_sorting_threshold=7)

    utils.assert_ordered_smiles_in_oeb_are_equal(sorted_by_fingerprint_oeb, \
            utils.get_list_of_canonical_isomeric_smiles(["N"] * 10 + ["C#N"] * 10 + ["O=C=O"] * 10))
    utils.assert_ordered_smiles_in_file_are_equal(smiles_dataset_file, \
            utils.get_list_of_canonical_isomeric_smiles(["N", "C#N", "O=C=O"]))
