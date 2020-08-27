"""
Script for filtering a dataset for molecules with single 
nitrogen-nitrogen bonds
Usage:
    python surrounding_n_n_pipeline.py --smiles-database database.smi
"""

import argparse
import n_n_fingerprint_funcs
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

    def relevant_if_has_single_bonded_nitrogen(mol: oechem.OEMol):
        return SUBS.SingleMatch(mol)

    def neighbors_and_wbo_fingerprint(mol: oechem.OEMol):
        match = SUBS.Match(mol, True)
        neighbors = n_n_fingerprint_funcs.find_neighboring_atoms(mol, match)
        wbo = n_n_fingerprint_funcs.wiberg_bond_order(mol, match)
        return neighbors + wbo

    config = {
        "database_type": "SMILES",
        "database_info": args.smiles_database,
        "relevance_function": relevant_if_has_single_bonded_nitrogen,
        "filter_output_oeb": "results/filter_output.oeb",
        "fingerprint_function": neighbors_and_wbo_fingerprint,
        "fingerprint_output_oeb": "results/fingerprint_output.oeb",
        "selection_frequency": 3,
        "dataset_type": "SMILES",
        "dataset_info": "results/dataset.smi",
        "sorted_by_fingerprint_oeb": "results/sorted_by_fingerprint.oeb",
        "in_memory_sorting_threshold": 25000,
        "tmpdir": "results/",
    }

    run_dance(config)


if __name__ == "__main__":
    main()
