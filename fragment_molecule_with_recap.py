#尝试使用RECAP进行分子碎裂！

#导入模块，载入分子
#!usr/bin/python3
from rdkit import Chem
from rdkit.Chem import Recap
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
mol = Chem.MolFromSmiles('O=c1nc([nH]c2NCC(Nc12)CN(c1ccc(cc1)C(=O)NC(CCC(=O)O)C(=O)O)C=O)N')
#结果以树结构形式获得
node = Recap.RecapDecompose(mol)
print (node)
#输出叶节点片段
leaves = [leaf.mol for leaf in node.GetLeaves().values()]
img = Draw.MolsToGridImage(leaves)
img.save('leaves.png')
#输出所有片段
all_nodes = [node.mol for node in node.GetAllChildren().values()]
img = Draw.MolsToGridImage(all_nodes)
img.save('all_nodes.png')