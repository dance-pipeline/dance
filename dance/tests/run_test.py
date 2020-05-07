"""Tests for run.py"""
import random

import pytest

from dance.run import run_dance

from . import utils

# pylint: disable=missing-function-docstring,invalid-name,unused-variable


def test_run_dance_raises_error_with_insufficient_config():
    with pytest.raises(RuntimeError):
        run_dance({})


def test_run_dance_correctly_completes_basic_pipeline(pipeline_test_files):
    """This is (mostly) the example described in the Getting Started documentation."""
    (smiles_file, filter_output_oeb, fingerprint_output_oeb, smiles_dataset_file,
     sorted_by_fingerprint_oeb) = pipeline_test_files

    smiles = ["N", "O=C=O", "C#N"] * 10
    random.shuffle(smiles)  # The pipeline should work for any molecule ordering.
    smiles_file.write_text("\n".join(smiles))

    def relevant_if_has_oxygen(mol):
        return any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())

    def size_fingerprint(mol):
        return (mol.NumAtoms(), )

    config = {
        "database_type": "SMILES",
        "database_info": smiles_file,
        "relevance_function": relevant_if_has_oxygen,
        "filter_output_oeb": filter_output_oeb,
        "fingerprint_function": size_fingerprint,
        "fingerprint_output_oeb": fingerprint_output_oeb,
        "selection_frequency": 3,
        "dataset_type": "SMILES",
        "dataset_info": smiles_dataset_file,
        "sorted_by_fingerprint_oeb": sorted_by_fingerprint_oeb,
        "in_memory_sorting_threshold": 6,
        "verbose": False,
    }

    run_dance(config)

    utils.assert_ordered_smiles_in_file_are_equal(smiles_dataset_file, \
            utils.get_list_of_canonical_isomeric_smiles(["O=C=O"] * 4))
