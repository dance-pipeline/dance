# Surrounding Nitrogen Molecules

This directory contains an example of using DANCE to select molecules with 
single nitrogen to nitrogen bonds. The main purpose of this directory is to
test a fingerprint that compares molecules based off the atoms surrounding 
their nitrogen-nitrogen bonds and the Wiberg bond orders of the bond(s). As
this fingeprint is still being tested and is computationally expensive, it is
not yet meant for larger databases like eMolecules.

## Instructions

### Main Pipeline

#### Dependencies

** You will need a cluster with slurm installed to run most of the steps below.**
Install the following dependencies:
- dance (see https://dance.readthedocs.io)
- openeye (`conda install -c openeye openeye-toolkits` or
  `pip install -i https://pypi.anaconda.org/openeye/simple openeye-toolkits`)
  - We used version 2019.10.2 in particular
- openforcefield (`conda install -c omnia openforcefield`)
  - We used version 0.6.0 in particular
  
Download the SMILES version of eMolecules from:
https://www.emolecules.com/info/plus/download-database

#### Preprocessing

By default, eMolecules contains salts along with the molecules. Run the
following command to generate `emolecules_no_salts.smi`, which has the salts
removed.

```bash
python remove_salts.py emolecules.smi emolecules_no_salts.smi
```

This command should take at most a few minutes to run.

As this file is not yet intended to run on the entirety of emolecules, it 
should be shortened to run quicker allowing for more efficient testing. Run 
the following terminal commands to only utilize the first 1000 molecules of 
the emolecules database.

```bash
tail -n+2 emolecules_no_salts.smi > emolecules_no_first_line.smi
```

```bash
head -n1000 emolecules_no_first_line.smi > nosalt_first1000.smi
```

Finally, create a local results directory and a tmp directory within it:

```bash
mkdir -p results/tmp
```

#### Running the Pipeline

Execute the following command to run the pipeline.

```bash
python surrounding_n_n_pipeline.py --smiles-database nosalt_first1000.smi
```

The final datasets will be found under the newly created results folder. This
command should take only a few seconds to run. To easily read the created OEB
files in a .txt format, run the following command:

```bash
 python oeb_results_to_txt.py results/
```

#### Results 

This pipeline begins with 1000 emolecules. The selection frequency of 3 
in the configuration for the pipeline allows for 334 potenital molecules to be
recorded in the results. However, after filtering through the first 1000
emolecules, the pipeline results in only26 molecules that were identified to 
have at least 1 single nitrogen-nitrogen bond. 

These 26 molecules are fingerprinted with by a list of 15 values. The first 
10 values in the list are reserved for atomic numbers of neighboring atoms to
the nitrogen-nitrogen bond(s) found in the given molecule. The remaining 5 
values in the list are reserved for Wiberg bond order values of the 
nitrogen-nitrogen bond(s) in the molecule. 

It is very likely a molecule's fingerprint data will not have a value for 
every index in the list. This is due to the possibility that molecules may not
have multiple nitrogen-nitrogen bonds. The unused values in a molecule's 
fingerprint are set to 0. 

In the instance a molecule is unable to be configured in order to calculate 
its Wiberg bond order, an error message will be printed including the smiles 
string for the molecule. The 5 values of the list set aside for Wiberg bond
order values will then be replaced with a -1. Currently, some molecules will
display a warning that they failed due to unspecified stereochemistry. This
message can be ignored and a molecule will be known to fail if the message,
"Failed to generate conformer for <molecule>" is displayed.