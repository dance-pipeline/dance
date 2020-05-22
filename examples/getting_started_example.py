"""Demo script that filters a SMILES database for molecules with oxygens.

Based on the example in Getting Started.

Usage:
    python getting_started_example.py --smiles-database database.smi
"""
import argparse

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

    def relevant_if_has_oxygen(mol):
        return any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())

    def size_fingerprint(mol):
        return (mol.NumAtoms(), )

    config = {
        "database_type": "SMILES",
        "database_info": args.smiles_database,
        "relevance_function": relevant_if_has_oxygen,
        "fingerprint_function": size_fingerprint,
        "selection_frequency": 3,
        "dataset_type": "SMILES",
        "dataset_info": "dataset.smi",
        "in_memory_sorting_threshold": 25000,
    }

    run_dance(config)


if __name__ == "__main__":
    main()
