from openforcefield.topology import Molecule
import qcengine
import qcelemental as qcel
from qcelemental.models import AtomicInput
from qcelemental.models.common_models import Model
from openeye import oequacpac, oechem

import numpy as np
import matplotlib.pyplot as plt

def get_xtb(mol):
    qc_mol = mol.to_qcschema()

    xtb_model = Model(method="gfn2-xtb", basis=None)
    qc_task = AtomicInput(molecule=qc_mol, driver="energy", model=xtb_model)

    result = qcengine.compute(input_data=qc_task, program="xtb")
    
    # xtb returns energy in hartree
    return result.return_result * qcel.constants.conversion_factor("hartree", "kcal/mol")

def get_am1(mol):
    oe_mol = mol.to_openeye()
    
    calc = oequacpac.OEAM1()
    result = oequacpac.OEAM1Results()
    calc.CalcAM1(result, oe_mol)
    # am1 returns energy in kcal/mol
    return result.GetEnergy()

def get_conformers(openff_mol, max_conformers):
    mol = openff_mol.to_openeye()
    
    x = []
    for i in mol.GetConfs():
        # convert each conformer back to openff molecule
        ii = oechem.OEMol(i)
        aa = Molecule.from_openeye(ii, allow_undefined_stereo=True)

        x.append(aa)
        if len(x) == max_conformers:
            break
    return x

def calc_rmsd(vals):
    return np.sqrt(np.nanmean(vals**2,axis=0))


if __name__ == "__main__":
    import sys

    if len(sys.argv) == 1:
        print("\tUsage: python3 xtb_am1_ratio.py data.smi [amount_of_molecules]")
        sys.exit(0)

    filename = sys.argv[1]
    
    number = float("inf")
    if len(sys.argv) > 2:
        try:
            number = int(sys.argv[2])
        except ValueError:
            print("\tUsage: python3 xtb_am1_ratio.py data.smi [amount_of_molecules]")
            print("\t\tamount_of_molecules must be int")
            sys.exit(1)
        
    if len(sys.argv) > 3:
        scale = bool(sys.argv[3])


    import os
    if not os.path.exists('plots'):
        os.makedirs('plots')


    print(f"Beginning benchmark on {filename}" + ("" if number == float("inf") else f", on the first {number} molecules"))

    print(f"Reading molecules")
    mols = []
    xtb_en = []
    am1_en = []
    ffmols = []

    cnt = 0
    skipped = 0
    pkldata = {}

    with open(filename) as inputs:
        for line in inputs:
            if cnt >= number:
                break

            sm_mol = line.split()[0]

            if sm_mol == "isosmiles":
                continue

            try:
                molecule = Molecule.from_smiles(sm_mol, allow_undefined_stereo=True)
                molecule.generate_conformers()
                ffmols.append(get_conformers(molecule, 10))
                mols.append(sm_mol)
                pkldata[sm_mol] = {}
            except:
                skipped += 1
                continue
            cnt += 1

            print(f"\rFinished: {cnt}", end="")

    mols = np.array(mols)

    print()
    print(f"Skipped: {skipped}")

    print(f"Computing energy ratios")
    b = np.full((10, len(ffmols)), np.nan)

    for i, confs in enumerate(ffmols):
        xtb_en = np.zeros((len(confs)))
        am1_en = np.zeros((len(confs)))
        for j, conf in enumerate(confs):
            xtb_en[j] = get_xtb(conf)
            am1_en[j] = get_am1(conf)

        ratio = xtb_en / am1_en

        xtb_rmsd = calc_rmsd(xtb_en)
        am1_rmsd = calc_rmsd(am1_en)
        xtb_mean = np.mean(xtb_en)
        am1_mean = np.mean(am1_en)
        pkldata[mols[i]] = {
                            "xtb_rmsd" : xtb_rmsd,
                            "am1_rmsd" : am1_rmsd,
                            "xtb_mean" : xtb_mean,
                            "am1_mean" : am1_mean,
                            "spread" : np.max(ratio) - np.min(ratio),
                            "num_conformers" : len(confs),
                            "conformer_xtb_energies" : xtb_en,
                            "conformer_am1_energies" : am1_en,
                            "conformer_energy_ratios" : ratio
                           }
        b[:len(confs),i] = ratio
        print(f"\rFinished: {i + 1}", end="")
    
    print()

    rmsd = calc_rmsd(b)
    rmsd_signed = rmsd.copy()
    for i in range(len(rmsd)):
        rmsd_signed[i] *= np.sign(pkldata[mols[i]]["xtb_mean"])*np.sign(pkldata[mols[i]]["am1_mean"])

        pkldata[mols[i]]["rmsd_energy_ratio"] = rmsd[i]
        pkldata[mols[i]]["rmsd_energy_ratio_signed"] = rmsd_signed[i]

    print("Generating pickle file")
    import pickle
    '''
    pickle format

    {
        smiles string of molecule :
        {
            "spread" : float,                       # the difference between the largest and smallest energy ratios 
            "num_conformers" : int,                 # the number of conformations sampled, up to 10
            "conformer_xtb_energies" : [float],     # array of xtb-calculated energies for each conformer
            "conformer_am1_energies" : [float],     # array of am1-calculated energies for each conformer
            "conformer_energy_ratios" : [float],    # energy ratios
            "rmsd_energy_ratio" : float             # the average energy ratio
        }
    }
    '''
    with open('xtb_am1_ratio.pickle', 'wb') as pkfile:
        pickle.dump(pkldata, pkfile)

    print("Generating plots")

    for i in range(int(np.ceil(len(mols)/20))):
        idxs = np.arange(i*20, min(len(mols), (i+1)*20) )
        for confidx in range(b.shape[0]):
            aa, = plt.plot(idxs, b[confidx,idxs], color="red", marker="_", markersize=15, linestyle="None")

        ab, = plt.plot(idxs,rmsd_signed[idxs], color="black", marker="x", markeredgewidth=2, linestyle="None")

        plt.legend((aa, ab), ("individual energy ratios", "signed RMSD energy ratio"))
        plt.title("Energy ratio am1 / xtb of up to 10 conformations of a molecule")
        plt.ylabel("energy ratio")
        plt.xticks(idxs)
        #plt.tight_layout()
        fig = plt.gcf()
        fig.set_size_inches(12.5, 6.5)
        plt.grid(axis="x")
        plt.savefig(f"plots/xtb_am1_ratio-{i}.png")#, bbox_inches="tight")
        print(f"Saving plots/xtb_am1_ratio-{i}.png")
        plt.close()
    

    with open("plots/molecule_legend.txt" ,'w') as lfile:
        for i, molname in enumerate(mols):
            lfile.write(f"{i}\t{molname}\n")

    print("Done")
