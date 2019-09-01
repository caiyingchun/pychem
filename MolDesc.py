import numpy as np
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

trainset = [mol for mol in Chem.SDMolSupplier('mols.sdf') if mol is not None]
nms=[x[0] for x in Descriptors._descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
trainDescrs = [calc.CalcDescriptors(x) for x in trainset]
trainDescrs = np.array(trainDescrs)
print trainDescrs
