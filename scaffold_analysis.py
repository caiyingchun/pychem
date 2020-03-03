#!/usr/bin/env python3

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold

drugbank_input = Chem.SDMolSupplier('drugbank.sdf')
drugbank = [m for m in drugbank_input if m]

basic_structure = drugbank[222]
atomic_scaffold = MurckoScaffold.GetScaffoldForMol(basic_structure)
atomic_scaffold.Compute2DCoords()
graph_scaffold = MurckoScaffold.MakeScaffoldGeneric(atomic_scaffold)
Draw.MolsToGridImage([basic_structure, atomic_scaffold, graph_scaffold])

drugbank_atomic_scaffolds = [MurckoScaffold.GetScaffoldForMol(mol) for mol in drugbank]
for i in drugbank_atomic_scaffolds:
    i.Compute2DCoords()
	
def genericize_scaffold(s):
    try:
        return MurckoScaffold.MakeScaffoldGeneric(s)
    except ValueError:
        return None
drugbank_grafh_scaffolds = [genericize_scaffold(s) for s in drugbank_atomic_scaffolds]

len(drugbank), len(drugbank_atomic_scaffolds), len(drugbank_grafh_scaffolds), len([x for x in drugbank_grafh_scaffolds if x == None])

Draw.MolsToGridImage([drugbank[111], drugbank_atomic_scaffolds[111], drugbank_grafh_scaffolds[111]])

scaffold_smiles = [Chem.MolToSmiles(scaffold) for scaffold in drugbank_grafh_scaffolds if scaffold != None]

len(scaffold_smiles), scaffold_smiles[111]

import collections
counter=collections.Counter(scaffold_smiles)
print(counter)

most_freq = Chem.MolFromSmiles('C1CCCCC1')
second_freq = Chem.MolFromSmiles('C1CCC(CC2CCCCC2)CC1')
Draw.MolsToGridImage([most_freq, second_freq])