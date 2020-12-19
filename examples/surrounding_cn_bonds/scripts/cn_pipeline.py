"""
Script for filtering a dataset for molecules with single 
nitrogen-nitrogen bonds
Usage:
    python cn_pipeline.py --smiles-database database.smi
"""

import os, sys

import argparse
from openeye import oechem

sys.path.append('/Users/williamyang/Documents/dance/dance')
from dance.run import run_dance

import cn_fingerprint_funcs

def main():
    """Builds and runs a pipeline."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles-database",
                        type=str,
                        required=True,
                        help=("Name of a SMILES file containing molecules "
                              "to process in the pipeline."))
    args = parser.parse_args()

    SUBS = oechem.OESubSearch("[#6:1]-[#7:2]")

    def relevant_if_has_single_cn_bond(mol: oechem.OEMol):
        match = SUBS.Match(mol, True)
        cnCount = 0

        while match.IsValid():
            for _ in match.Target().GetTargetBonds():
                cnCount += 1
            match.Next()

        if cnCount == 1:
            mol, status = cn_fingerprint_funcs.smiles_to_oemol(oechem.OEMolToSmiles(mol))
            return status
        
        return False

    def neighbors_and_wbo_fingerprint(mol: oechem.OEMol):
        result = []
        mol, status = cn_fingerprint_funcs.smiles_to_oemol(oechem.OEMolToSmiles(mol))
        
        if status is False:
            return [-1,-1,-1,-1,-1]
        
        match = SUBS.Match(mol, True)
        
        wbo = cn_fingerprint_funcs.wiberg_bond_order(mol, match)
        result.append(wbo)
        
        neighbors = cn_fingerprint_funcs.find_neighboring_atoms(mol, match)
        result += neighbors
    def size_fingerprint(mol):
        return (mol.NumAtoms(),)
        
        
        return result

    config = {
        "database_type": "SMILES",
        "database_info": args.smiles_database,
        "relevance_function": relevant_if_has_single_cn_bond,
        "filter_output_oeb": "results/filter_output.oeb",
        "fingerprint_function": size_fingerprint,
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