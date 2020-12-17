"""
Script for filtering a dataset for molecules with single 
nitrogen-nitrogen bonds
Usage:
    python surrounding_nn_pipeline.py --smiles-database database.smi
"""

import argparse
import nn_fingerprint_funcs
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
        match = SUBS.Match(mol, True)
        nnCount = 0

        while match.IsValid():
            for i in match.Target().GetTargetBonds():
                nnCount += 1
            match.Next()

        if nnCount == 1:
            mol, status = nn_fingerprint_funcs.smiles2oemol(oechem.OEMolToSmiles(mol))
            return status
        
        return False

    def neighbors_and_wbo_fingerprint(mol: oechem.OEMol):
        result = []
        mol, status = nn_fingerprint_funcs.smiles2oemol(oechem.OEMolToSmiles(mol))
        
        if not status:
            return [-1,-1,-1,-1,-1]
        
        match = SUBS.Match(mol, True)
        
        wbo = nn_fingerprint_funcs.wiberg_bond_order(mol, match)
        result.append(wbo)
        
        neighbors = nn_fingerprint_funcs.find_neighboring_atoms(mol, match)
        result += neighbors
        
        
        return result

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
