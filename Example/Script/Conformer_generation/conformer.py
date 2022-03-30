## Import module
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign


################################################################
## This is used to analyze docking result 
## The result is saved as a csv file. 
## Plot :	top rank ranking accuracy	
##		distribution of # of rotors
################################################################


################################################################
## Bash -- prepare folder
################################################################

subprocess.call(['rm', '-r', '-f', 'conformer'])
subprocess.call(['mkdir', '-p', 'conformer'])

################################################################
## Read in ligand and generate conformers 
################################################################

## Read in ligand
mol = Chem.MolFromMolFile('ligand.mol')
mol = Chem.AddHs(mol)

## Generate conformers
params = AllChem.ETKDGv3()
params.useRandomCoords = True
params.useSmallRingTorsions = True
params.useMacrocycleTorsions = True
params.useExpTorsionAnglePres = True
cutoff = np.linspace(2, 0.5, num = 16)
for radius in cutoff :
	params.pruneRmsThresh = radius
	cids = AllChem.EmbedMultipleConfs(mol, numConfs = changeNum, params = params)
	print("With raidus cutoff " + str(radius) + ", we have " + str(len(cids)) + " conformers generated")
	if len(cids) > 10 : break
rdMolAlign.AlignMolConformers(mol)

## Save conformers
for cid in cids :
	w = Chem.PDBWriter('tmp.pdb')
	w.write(mol, confId = cid)
	w.close()
	idx = str(cid + 1)
	subprocess.call(['bash', 'rename.sh', idx])

exit()

#!-------------------------stop here-----------------------------

#for cid in cids:
#	w = Chem.PDBWriter('tmp.pdb')
#	w.write(mol, confId = cid)
#	w.close
#	idx = str(cid + 1)
#	subprocess.call(['bash', 'save.sh', idx])

#res = []
#mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant = 'MMFF94s')
#for cid in cids :
#	ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId = cid)
#	e = ff.CalcEnergy()
#	res.append((cid, e))

#sorted_res = sorted(res, key = lambda x:x[1])

#for cid, e in sorted_res :
for cid in cids :
	w = Chem.PDBWriter('tmp.pdb')
	#mol.SetProp('CID', str(cid))
	#mol.SetProp('Energy', str(e))
	w.write(mol, confId = cid)
	w.close()
	idx = str(cid + 1)
	subprocess.call(['bash', 'rename.sh', idx])

