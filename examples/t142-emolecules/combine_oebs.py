#!/usr/bin/env python3

#SBATCH --job-name=save_molecule_parameters
#SBATCH --partition=ilg2.3  # Modify as needed
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=30:00:00
#SBATCH --distribution=block:cyclic
#----------------------------------------
"""Combines multiple OEB files together into a final file.

The files are the ones output by running the `save_molecule_parameters.py` job
with slurm.

Usage:
    sbatch combine_oebs.py [NUM_FILES]
"""
import logging
import os
import sys

from openeye import oechem


def main():
    """Combines the files together."""
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

    if len(sys.argv) != 2:
        print("Usage: python combine_oebs.py [NUM_FILES]")
        sys.exit(1)

    num_files = int(sys.argv[1])
    output_file = "results/combined_molecules_with_parameters.oeb"
    ofs = oechem.oemolostream(output_file)

    total_molecules = 0

    for i in range(num_files):
        filename = f"results/molecules_with_params_{i}.oeb"
        logging.info("Reading molecules from %s", filename)
        ifs = oechem.oemolistream(filename)
        for mol in ifs.GetOEMols():
            oechem.OEWriteMolecule(ofs, mol)
            total_molecules += 1
        ifs.close()
    logging.info("Wrote %d molecules to %s", total_molecules, output_file)
    ofs.close()


if __name__ == "__main__":
    main()
