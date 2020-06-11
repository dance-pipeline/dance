#!/usr/bin/env python3

#SBATCH --job-name=save_molecule_parameters
#SBATCH --partition=ilg2.3  # Modify as needed
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=30:00:00
#SBATCH --distribution=block:cyclic
#----------------------------------------
"""Script for selecting a dataset with t142 params from a molecule database.

Usage:
    sbatch smirnoff_pipeline.py

Modify parameters below as needed.

    isort:skip_file
"""

import logging
import os
import sys

from openeye import oechem

# Necessary to add cwd to path when script run by slurm (since it executes
# a copy)
sys.path.append(os.getcwd())

# pylint: disable=wrong-import-position
import smirnoff_param_utils
from dance.run import run_dance


def main():
    """Builds and runs a pipeline."""
    logging.getLogger("dance").setLevel(logging.DEBUG)
    logging.getLogger("dance").addHandler(logging.StreamHandler(sys.stderr))

    def relevant_if_has_t142(mol: oechem.OEMol):
        params = smirnoff_param_utils.read_params_from_mol(mol)
        return "t142" in params

    def size_and_wbo_fingerprint(mol: oechem.OEMol):
        """Fingerprint with number of atoms in the mol and WBO between central atoms in the t142 param.

        (WBO is Wiberg Bond Order.)

        Note: The indices are the same in OpenEye and OFF molecules, so the
        parameter indices are correct, even though they were calculated for the
        OFF mol.
        """
        oechem.OEAddExplicitHydrogens(mol)
        params = smirnoff_param_utils.read_params_from_mol(mol)
        return (
            mol.NumAtoms(),
            smirnoff_param_utils.calculate_t142_central_wbo(mol, params),
        )

    config = {
        "database_type": "OEB",
        "database_info": "combined_molecules_with_parameters.oeb",
        "relevance_function": relevant_if_has_t142,
        "filter_output_oeb": "results/filter_output.oeb",
        "fingerprint_function": size_and_wbo_fingerprint,
        "fingerprint_output_oeb": "results/fingerprint_output.oeb",
        "selection_frequency": 375,  # To select ~20 molecules.
        "dataset_type": "OEB",
        "dataset_info": "results/dataset.oeb",
        "sorted_by_fingerprint_oeb": "results/sorted_by_fingerprint.oeb",
        "in_memory_sorting_threshold": 25000,
        "tmpdir": "results/tmp",
    }

    run_dance(config)


if __name__ == "__main__":
    main()
