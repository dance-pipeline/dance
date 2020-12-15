'''
Script to cleanly portray the fingerprint data of molecules

Usage:
    python get_fingerprint_data.py
'''

from openeye import oechem

def read():
    ifs = oechem.oemolistream("results/sorted_by_fingerprint.oeb")
    ifs.SetFormat(oechem.OEFormat_OEB)
    
    molList = []
    molCount = 1
    
    for mol in ifs.GetOEGraphMols():
        molList.append(oechem.OEGraphMol(mol))
        
    with open("results/sorted_fingerprint_data.txt", "w") as data:
        with open("results/failed_molecules.txt", "w") as failedMols:
            with open("results/sorted_by_fingerprint.txt") as mols:
                for lineNum, line in enumerate(mols):
                    data.write(f"Mol #{molCount} ")
                    data.write(f"{line}")
                    dataVal = molList[lineNum].GetData()
                    if molList[lineNum].GetData("dance_fingerprint_value_4") == -1:
                        failedMols.write(f"{line}")
                    data.write(f"Fingerprint data: {dataVal}\n")
                    data.write("\n")
                    molCount += 1
                

if __name__ == '__main__':
    read()

