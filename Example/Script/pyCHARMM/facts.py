################################################################
## 
## pyCHARMM covalent docking facts implicit solvent minimization 
## Created by Yujin Wu, wyujin@umich.edu
## Date: 03/27/2022
##
################################################################

## Import module
import pycharmm
import pycharmm.psf as psf
import pycharmm.read as read
import pycharmm.lingo as lingo
import pycharmm.generate as gen
import pycharmm.energy as energy
import pycharmm.cons_fix as cons_fix
import pycharmm.minimize as minimize
import pycharmm.settings as settings
import pycharmm.implicit_solvent as implicit
import numpy as np
import pandas as pd
from os import listdir

## Topology and parameter files 
settings.set_bomb_level(-1)
read.rtf('"../Toppar/top_all36_prot.rtf"')
read.rtf('"../Toppar/top_all36_cgenff.rtf"', append = True)
read.prm('"../Toppar/par_all36m_prot.prm"', flex = True)
read.prm('"../Toppar/par_all36_cgenff.prm"', append = True, flex = True)
settings.set_bomb_level(0)
lingo.charmm_script('stream "./ligandrtf"')

## Build system
read.sequence_pdb("ligand.pdb")
settings.set_bomb_level(-1)
gen.new_segment(seg_name = "LIGA")
settings.set_bomb_level(0)
read.pdb("ligand.pdb", resid = True)
read.psf_card("protein.psf", append = True)
read.pdb("protein.pdb", resid = True)

## Set parameter
lingo.charmm_script('''faster on''')
nbond = pycharmm.NonBondedScript(
    nbxmod = 5, atom = True, cdiel = True, eps = 1,
    shift = True, vatom = True, vdistance = True, vswitch = True,
    cutnb = 14.0, ctofnb = 12.0, ctonnb = 10.0, e14fac = 1.0, wmin = 1.5).run()
lingo.charmm_script('''scalar wmain = radius''')
facts = implicit.FACTS(tcps = 22, teps = 1, gamm = 0.015, 
    tavw = True, conc = 0.1, temp = 298).run()
    
## Cons fix atoms
ligand_atom = np.loadtxt("ligandatom", dtype = str)
ligand = pycharmm.SelectAtoms().by_seg_id("LIGA")
for idx in np.arange(len(ligand_atom)) :
    if idx == 0 :
        ligand_warhead = pycharmm.SelectAtoms().by_atom_type(ligand_atom[idx])
    else :
        ligand_warhead = ligand_warhead | pycharmm.SelectAtoms().by_atom_type(ligand_atom[idx])
ligand_warhead = ligand & ligand_warhead 
receptor = ligand.__invert__()
fix_atoms = receptor | ligand_warhead 
print("Cons fix atoms")
print("Num of ligands", len(listdir(dock_result/run/cluster)))

## Loop through docked pose to get system FACTS rescored energy
total = []
ligandID = []
for idx in np.arange(len(listdir("dock_result/run/cluster"))) + 1 :
    read.pdb("dock_result/run/cluster/top_" + str(idx) + ".pdb", resid = True)
    cons_fix.setup(fix_atoms)

    ## FACTS MD
    minimize.run_abnr(nstep = 1000, tolgrd = 0.001)
    cons_fix.turn_off()

    ## Get system energy
    energy.show()
    total.append(energy.get_total())
    ligandID.append("top_" + str(idx))

## Set parameter
psf.delete_atoms(selection = receptor)
lingo.charmm_script('''faster on''')
nbond = pycharmm.NonBondedScript(
    nbxmod = 5, atom = True, cdiel = True, eps = 1,
    shift = True, vatom = True, vdistance = True, vswitch = True,
    cutnb = 14.0, ctofnb = 12.0, ctonnb = 10.0, e14fac = 1.0, wmin = 1.5).run()
lingo.charmm_script('''scalar wmain = radius''')

## Loop through docked pose to get ligand FACTS rescored energy
ligandEner = []
for idx in np.arange(len(listdir("dock_result/run/cluster"))) + 1 :
    read.pdb("dock_result/run/cluster/top_" + str(idx) + ".pdb", resid = True)
    minimize.run_abnr(nstep = 1000, tolgrd = 0.001)
    ligandEner.append(energy.get_total())

## Get covalent bond energy
dockresult = pd.read_csv("dock_result/run/clusterResult.csv", sep = "	")
covBond = dockresult["grid_hbond"].to_numpy()

## Summary result
total = np.asarray(total) + covBond
ligandEner = np.asarray(ligandEner)
facts_result = pd.DataFrame()
facts_result['Total_energy'] = total
facts_result['Covalent_bond'] = covBond
facts_result['Ligand_energy'] = ligandEner
facts_result['pdbID'] = ligandID
facts_result.to_csv('dock_result/run/factsResult.csv', sep = '	')

