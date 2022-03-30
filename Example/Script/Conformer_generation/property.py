## Import module
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np


################################################################
## This is used to analyze ligand property 
################################################################


################################################################
## Read in ligand 
################################################################

## Read in ligand
mol = Chem.MolFromMolFile('ligand.mol')
mol = Chem.AddHs(mol)

## Find rotatable bond 
'''
	This is used to find all rotatable bond
	Get rid of terminal rotatable bonds
'''

rotatableBond = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
match = mol.GetSubstructMatches(rotatableBond)
total = np.shape(match)[0]
length = np.linspace(0, total - 1, total).astype(int)
array = []
count = 1

for idxI in length :
	terminal = False

	atom = mol.GetAtomWithIdx(match[idxI][0])
	neighbor = np.array([x.GetAtomicNum() for x in atom.GetNeighbors()])
	if (neighbor == 1).sum() == np.shape(neighbor)[0] - 1 : terminal = True

	atom = mol.GetAtomWithIdx(match[idxI][1])
	neighbor = np.array([x.GetAtomicNum() for x in atom.GetNeighbors()])
	if (neighbor == 1).sum() == np.shape(neighbor)[0] - 1 : terminal = True

	if not terminal :
		tmp = np.array([match[idxI][0], match[idxI][1]])
		if count == 1 : 
			array = tmp 
		else : 
			array = np.concatenate((array, tmp), axis = 0)
		count += 1


## Classify rotatable bond
'''
	This is used to classify ligand rotatable bond
	total 	--> the total number of the rotatable bonds after get rid of terminal rotatable bond
	free	--> count of freely rotatable bond
	conj	--> count of conjugated rotatable bond
'''

free = 0
if np.ndim(array) == 1 : array = np.reshape(array, (-1, 2))
total = np.shape(array)[0]
length = np.linspace(0, total - 1, total).astype(int)
for idxI in length :
	if (array == array[idxI, 0]).sum() == 1 and (array == array[idxI, 1]).sum() == 1 : free += 1
conj = total - free

## Output result
print(free, conj, Descriptors.MolWt(mol))

