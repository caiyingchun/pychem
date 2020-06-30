import sys
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

#SMILES = 'O=C(NCCCCCCCCNC(C[C@H]1CNC=CC=N1)=O)COC2=CC=CC=C2'

HowMany = 5000
limitRMS = 0.2

#mol = Chem.MolFromSmiles(SMILES)
mol = Chem.SDMolSupplier('linker.sdf')[0]
mol = AllChem.AddHs(mol)
allconf = AllChem.EmbedMultipleConfs(mol,
   numConfs = HowMany, enforceChirality=True, pruneRmsThresh = limitRMS)
w = Chem.SDWriter('confs_rmsd0.2.sdf')
for confId in allconf:
   print(confId)
   ff = AllChem.UFFGetMoleculeForceField(mol, confId = confId)
   ff.Minimize()
   energy_value = ff.CalcEnergy()
   mol.SetProp('_Name', 'linker')
   mol.SetProp('ENERGY', '{0:.2f}'.format(energy_value))
   w.write(mol, confId = confId)
w.close()
