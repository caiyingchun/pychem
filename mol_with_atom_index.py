import copy
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

def mol_with_atom_index1(mol):
    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))
    return mol
    
def mol_with_atom_index2(mol):
    for atom in mol.GetAtoms():
        atom.SetProp('atomLabel',str(atom.GetIdx()))
    return mol