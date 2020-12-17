"""
Script for generating heatmaps based on the MACCS key fingerprint values of
nitrogen-nitrogen molecules
"""

import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import nn_fingerprint_funcs
import numpy as np
from openeye import oechem, oegraphsim

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles_database",
                        type=str,
                        required=True,
                        help=("Name of a SMILES file containing molecules "
                              "to generate a MACCS key fingerprint heatmap"))
    args = parser.parse_args()
    
    ifs = oechem.oemolistream(args.smiles_database)
    
    molecules = get_grouped_mols(ifs)
    
    for mol_group in molecules:
        maccs_keys = get_MACCS_keys(mol_group)
        generate_heatmap(maccs_keys, mol_group, molecules.index(mol_group)+1)

def get_MACCS_keys(molecules):
    """
    Obtains and returns a tanimoto comparison of the MACCS key fingeprint values for a list of molecules
    """
    maccs_keys = []
    for molA in molecules:
        tanimoto_vals = []
        fpA = oegraphsim.OEFingerPrint()
        oegraphsim.OEMakeFP(fpA, molA, oegraphsim.OEFPType_MACCS166)
        
        # Compares molA against every other mol in the list
        for molB in molecules:
            fpB = oegraphsim.OEFingerPrint()
            oegraphsim.OEMakeFP(fpB, molB, oegraphsim.OEFPType_MACCS166)
            val = oegraphsim.OETanimoto(fpA, fpB)
            tanimoto_vals.append( round(val, 2) )
        maccs_keys.append(tanimoto_vals)
    
    return maccs_keys
    
def get_grouped_mols(ifs):
    """
    Sorts the molecules in the provided smiles_database based on number of atoms and groups the molecules into groups with a max of 25 per group for the purpose of creating better looking plots
    """
    molecules = []
    
    for molecule in ifs.GetOEMols():
        mol, status = nn_fingerprint_funcs.smiles2oemol(oechem.OEMolToSmiles(molecule))
        if (status):
            molecules.append(mol)
    
    molecules.sort(key=lambda mol: mol.NumAtoms())
    
    groups_of_25 = [molecules[i:i+25] for i in range(0, len(molecules), 25)]
        
    return groups_of_25
    
    
def generate_heatmap(maccs_keys, mols, group_num):
    """
    Generates a heatmap using the tanimoto values collected by comparing the MACCS key fingerprint values of molecules
    """
    map = np.array(maccs_keys)

    fig, ax = plt.subplots(figsize = (10,10))
    im = ax.imshow(map, cmap="YlGn")

    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_ylabel(ylabel = "", rotation=-90, va="bottom")

    ax.set_xticks(np.arange(len(mols)))
    ax.set_yticks(np.arange(len(mols)))

    ax.set_xticklabels([num for num in range(len(mols))])
    ax.set_yticklabels([num for num in range(len(mols))])

    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
                       
    plt.setp(ax.get_xticklabels(), rotation=0, ha="right", rotation_mode="anchor")

    for i in range(len(mols)):
        for j in range(len(mols)):
            text = ax.text(j, i, maccs_keys[i][j],
                           ha="center", va="center", color="black")
                           
    fig.tight_layout()
    plt.title(f"N-N Molecule MACCS Key Fingerprint Heatmap (Mol group {group_num}: {len(mols)} mols)")
    # plt.show()
    fig.savefig(f"analysis/MolGroup{group_num}.png", dpi = fig.dpi)

if __name__ == "__main__":
    main()
