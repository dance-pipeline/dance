"""
Script for filtering a dataset for molecules with single 
nitrogen-nitrogen bonds

Usage:
    python n_n_pipeline.py --smiles-database database.smi
"""

import argparse
from openeye import oechem
from dance.run import run_dance


def main():
    """Builds and runs a pipeline."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles-database",
                        type=str,
                        required=True,
                        help=("Name of a SMILES file containing molecules "
                              "to process in the pipeline."))
    args = parser.parse_args()

    SUBS = oechem.OESubSearch("[#7:1]-[#7:2]")

    def relevant_if_has_single_bonded_nitrogen(mol):
        return SUBS.SingleMatch(mol)

    def size_fingerprint(mol):
        return (mol.NumAtoms(), )

    config = {
        "database_type": "SMILES",
        "database_info": args.smiles_database,
        "relevance_function": relevant_if_has_single_bonded_nitrogen,
        "fingerprint_function": size_fingerprint,
        "selection_frequency": 3,
        "dataset_type": "SMILES",
        "dataset_info": "dataset.smi",
        "in_memory_sorting_threshold": 25000,
    }

    run_dance(config)


if __name__ == "__main__":
    main()
