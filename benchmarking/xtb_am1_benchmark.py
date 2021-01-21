from openforcefield.topology import Molecule
import qcengine
import qcelemental as qcel
from qcelemental.models import AtomicInput
from qcelemental.models.common_models import Model
from openeye import oequacpac, oechem
import matplotlib.pyplot as plt



def get_xtb(mol):
    qc_mol = mol.to_qcschema()

    xtb_model = Model(method="gfn2-xtb", basis=None)
    qc_task = AtomicInput(molecule=qc_mol, driver="energy", model=xtb_model)

    result = qcengine.compute(input_data=qc_task, program="xtb")

    # xtb returns energy in hartree
    return result.return_result * qcel.constants.conversion_factor("hartree", "kcal/mol")

import math

def get_am1(mol):
    oe_mol = mol.to_openeye()
    calc = oequacpac.OEAM1()
    result = oequacpac.OEAM1Results()
    calc.CalcAM1(result, oe_mol)
    # am1 returns energy in kcal/mol
    return result.GetEnergy()

def getEnergies(sm_mol : str) -> ("xtb energy", "am1 energy"):
    molecule = 0

    try:
        molecule = Molecule.from_smiles(sm_mol, allow_undefined_stereo=True)
        molecule.generate_conformers()
    except Exception:
        return (None, None)

    return (get_xtb(molecule), get_am1(molecule))



if __name__ == "__main__":
    import sys

    if len(sys.argv) == 1:
        print("\tUsage: python3 xtb_am1_benchmark.py data.smi [amount_of_molecules]")
        sys.exit(0)

    filename = sys.argv[1]
    
    number = float("inf")
    if len(sys.argv) > 2:
        try:
            number = int(sys.argv[2])
        except ValueError:
            print("\tUsage: python3 xtb_am1_benchmark.py data.smi [amount_of_molecules]")
            print("\t\tamount_of_molecules must be int")
            sys.exit(1)
        

    mols = []
    xtb_en = []
    am1_en = []
    data = {}

    cnt = 0

    print(f"Beginning benchmark on {filename}" + ("" if number == float("inf") else f", on the first {number} molecules"))

    with open(filename) as inputs:
        for line in inputs:
            if cnt >= number:
                break

            sm_mol, v_id, p_id = line.split()

            if sm_mol == "isosmiles":
                continue

            mols.append(sm_mol)
            _xtb, _am1 = getEnergies(sm_mol)
            xtb_en.append(_xtb)
            am1_en.append(_am1)
            data[sm_mol] = {"xtb": _xtb, "am1": _am1}

            cnt += 1

            print(f"\rFinished: {cnt}", end="")

    print()
    print("Generating pickle file")

    import pickle
    with open('xtb_am1_benchmark.pickle', 'wb') as pkfile:
        pickle.dump(data, pkfile)

    print("Pickle file generated at xtb_am1_benchmark.pickle")
        
    print("Plotting data")

    scat_xtb = plt.scatter(mols, xtb_en)
    scat_am1 = plt.scatter(mols, am1_en)
    plt.legend((scat_xtb, scat_am1), ("xtb", "am1"))
    plt.ylabel("energy, kcal/mol")
    plt.xticks(rotation=90)
    fig = plt.gcf()
    #fig.set_size_inches(18.5, 10.5)
    plt.grid(axis="x")
    plt.show()

    print("Done")
