"""Pytest directory configuration file.

Fixtures defined here will be available to all tests in the directory.
"""
import pytest


@pytest.fixture
def pipeline_test_files(tmp_path):
    """Fixture for easily creating files for testing.

    Creates a SMILES file, an OEB filepath for the filter step, an OEB filepath
    for the fingerprint step, a SMILES filepath for outputting the final
    dataset, and an OEB filepath for storing the molecules sorted by
    fingerprint.
    """
    smiles_file = tmp_path / "smiles.smi"
    filter_output_oeb = tmp_path / "filter_output.oeb"
    fingerprint_output_oeb = tmp_path / "fingerprint_output.oeb"
    smiles_dataset_file = tmp_path / "smiles_output.smi"
    sorted_by_fingerprint_oeb = tmp_path / "sorted_by_fingerprint.oeb"
    return smiles_file, filter_output_oeb, fingerprint_output_oeb, smiles_dataset_file, sorted_by_fingerprint_oeb
