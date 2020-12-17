'''
Script to read .oeb files and write them to a readable .txt format

Usage:
    python oeb_results_to_txt.py results/
'''

from openeye import oechem
import glob
import sys

def read():
    for file in glob.glob(f"{sys.argv[1]}*.oeb"):
        
        ifs = oechem.oemolistream(file)
        ifs.SetFormat(oechem.OEFormat_OEB)
        ofs = oechem.oemolostream(f'{file[:-4]}.txt')
    
        for mol in ifs.GetOEGraphMols():
            oechem.OEWriteMolecule(ofs, mol)

if __name__ == '__main__':
    read()