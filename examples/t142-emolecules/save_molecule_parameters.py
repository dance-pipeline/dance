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
"""Slurm script for running parameter calculations.

Usage:
    Set the parameters in this script (see below), and run with:

        sbatch --array=0-{NUM_JOBS-1} save_molecule_parameters.py

    The script will select a random subset of the molecules in the input file,
    and each of the tasks in the job array will calculate parameters for a
    portion of the subset. In the command above, {NUM_JOBS-1} is the number of
    jobs to split the calculations into. Make sure to set NUM_JOBS in the
    parameters below as well.

    WARNING: Order matters in the command above. Do not place the `--array`
    after `save_molecule_parameters.py`.

    isort:skip_file
"""

import logging
import os
import random
import sys

from openeye import oechem

# Necessary to add cwd to path when script run by slurm (since it executes
# a copy)
sys.path.append(os.getcwd())

# pylint: disable=unused-import, wrong-import-position
import smirnoff_param_utils
from smirnoff_param_utils import \
    count_total_molecules

# ========================== PARAMETERS ==========================

# Name of the input file with the molecule database.
INPUT_FILE = "emolecules_no_salts.smi"

# Number in eMolecules -- if unknown, set to `count_total_molecules(INPUT_FILE)`
TOTAL_MOLECULES = 25946988

# Number of molecules to calculate parameters for.
NUM_MOLS_TO_CALCULATE = 10**6

# Number of jobs across which to split the calculations.
NUM_JOBS = 20

# Seed for generating indices -- currently set to the meaning of life, the
# universe, and everything ;)
RANDOM_SEED = 42

# ================================================================


def read_index_mols_from_file(filename: str, indices: set) -> oechem.OEMol:
    """Generates the molecules at the given indices in the given file"""
    ifs = oechem.oemolistream(filename)
    index = 0
    for mol in ifs.GetOEMols():
        if index in indices:
            yield oechem.OEMol(mol)
        index += 1


def main():
    """Calculates parameters on desired molecules and saves them to the output file."""
    logging.getLogger().setLevel(logging.INFO)
    logging.info("===== Parameter Calculation SLURM JOB =====")
    logging.info("")
    logging.info("The job will be started on the following node(s):")
    logging.info(os.environ["SLURM_JOB_NODELIST"])
    logging.info("")
    logging.info("Slurm User:         %s", os.environ['SLURM_JOB_USER'])
    logging.info("Run Directory:      %s", os.getcwd())
    logging.info("Job ID:             %s", os.environ['SLURM_JOB_ID'])
    logging.info("Job Name:           %s", os.environ['SLURM_JOB_NAME'])
    logging.info("Partition:          %s", os.environ['SLURM_JOB_PARTITION'])
    logging.info("Number of nodes:    %s", os.environ['SLURM_JOB_NUM_NODES'])
    logging.info("Number of tasks:    %s", os.environ['SLURM_NTASKS'])
    logging.info("Submitted From:     %s", os.environ['SLURM_SUBMIT_HOST'])
    logging.info("Submit directory:   %s", os.environ['SLURM_SUBMIT_DIR'])
    logging.info("===== Parameter Calculation SLURM JOB =====")
    logging.info("")

    # Since we are running an array job, each job is assigned a task id.
    task_id = int(os.environ["SLURM_ARRAY_TASK_ID"])

    # Indices of molecules to take out of all the molecules in the input file.
    indices = random.sample(range(TOTAL_MOLECULES), \
        NUM_MOLS_TO_CALCULATE)[NUM_MOLS_TO_CALCULATE // NUM_JOBS * task_id : \
                               NUM_MOLS_TO_CALCULATE // NUM_JOBS * (task_id + 1)]
    logging.info("Indices: %s", indices)

    # File for storing the final molecules.
    output_oeb = f"results/molecules_with_params_{task_id}.oeb"

    output_stream = oechem.oemolostream(output_oeb)
    mols_processed = 0
    for mol in read_index_mols_from_file(INPUT_FILE, indices):
        if mols_processed % 10 == 0:
            logging.info("Processing molecule %d", mols_processed)
        mols_processed += 1

        params = smirnoff_param_utils.calculate_mol_params(mol)
        smirnoff_param_utils.write_params_to_mol(mol, params)
        oechem.OEWriteMolecule(output_stream, mol)


if __name__ == "__main__":
    main()
