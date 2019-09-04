from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

def MorganFP(smi, radius=3, nBits=2048):
    mol = Chem.MolFromSmiles(smi)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nBits)
    return [bit for bit in fp]

def AtomPairFP(smi, maxLength=30, nBits=2048):
    mol = Chem.MolFromSmiles(smi)
    fp = AllChem.GetHashedAtomPairFingerprintAsBitVect(mol, maxLength=maxLength, nBits=nBits)
    return [bit for bit in fp]

def TTFP(smi, nBits=2048):
    mol = Chem.MolFromSmiles(smi)
    fp = AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=nBits)
    return [bit for bit in fp]

def MolDesc(smi):
    mol = Chem.MolFromSmiles(smi)
    nms=[x[0] for x in Descriptors._descList]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
    return list(calc.CalcDescriptors(mol))

def fpDensity(fp):
    l = len(fp)
    return (l - fp.count(0)) / float(l)

def test():
    smi = 'O=C(OC)C1=CC=CC=C1'
    fp = MorganFP(smi)
    print len(fp)

if __name__ == '__main__':
    test()